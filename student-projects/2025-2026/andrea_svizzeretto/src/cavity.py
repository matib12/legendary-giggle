''' ### CAVITY ENVIRONMENT ### '''

import math
from typing import Optional, Tuple, Union
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque

import numpy as np
from scipy import constants as const
import oreonspy.utils as ut

import gymnasium as gym
from gymnasium import logger, spaces
from gymnasium.envs.classic_control import utils
from gymnasium.error import DependencyNotInstalled
from env_lib.rewards import get_reward_function
from env_lib.terminations import get_termination_condition

plt.rcParams['figure.figsize'] = [10, 6]

import oreonspy as op


class CavityEnv(gym.Env[np.ndarray, Union[int, np.ndarray]]):
    
    metadata = {
        "render_modes": ["human", "rgb_array"],
        "render_fps": 50,
    }

    def resize_action_space(self, action):
        if action > self.action_box_limit:
            self.action_space = spaces.Box(-np.abs(action), np.abs(action) + 1e-10, shape=(1,), dtype=np.float32)
        elif action < -self.action_box_limit:
            self.action_space = spaces.Box(-np.abs(action)-1e-10, np.abs(action), shape=(1,), dtype=np.float32)

    def v_cr(self, cavity, lambd):
        return lambd / (2. * cavity.Finesse() * cavity.tau())  # [m/s]
    
    def time_window(self, v, cavity, lambd, number_of_FSR=1.):
        t_stop = number_of_FSR * lambd/(2.*v)
        number_of_points = int(np.ceil(t_stop/cavity.Theta))
        return number_of_points, np.linspace(0., t_stop, number_of_points)
    
    
    def set_eval_mode(self, is_eval=True):
        ''' Set the environment to evaluation mode.
         You have to call this method after unwrapping the environment and before starting 
         the evaluation episodes with evaluate_policy function.'''
        self.eval_mode = is_eval

    def __init__(self, render_mode: Optional[str] = None, reward_fn = "Modified_NewReward_ActionAware", term_fn = "simple_termination", cavity_to_load= None, 
                 
                 t_a = 0.1, r_a = 0.9, r_b = 0.9, L=3000., f_c = 150.e3, E_in = 1, lambd = 1064e-9, v_max = 1, random_L = False, max_steps = None,
                 
                 displacement_factor = 1., in_resonance_L = False, no_render=False, P_noise_sigma=0., PDH_noise_sigma=0.):

        self.no_render = no_render
        self.cavity_to_load = cavity_to_load
        # LASER
        self.E_in = E_in  
        self.lambd = lambd  # [m]
        self.k = 2.*np.pi / lambd

        # CAVITY
        self.__L__ = 0.
        self.in_resonance_L = in_resonance_L
        if self.in_resonance_L:
            self.__L__ = np.ceil(L/self.lambd)*self.lambd
        else:
            self.__L__ = L
            
        self.cavity = op.Cavity(T_a = t_a , R_a =r_a, R_b=r_b , L=self.__L__)

        if self.cavity_to_load is not None:   # If a cavity is provided, load it and set its length to the desired one (in resonance or not)
            self.cavity = self.cavity_to_load
            if self.in_resonance_L:
                self.__L__ = np.ceil(self.cavity.__L__/self.lambd)*self.lambd
                self.cavity.__L__ = self.__L__
            else:
                self.__L__ = self.cavity.__L__

        # SAMPLING 
        self.f_c = f_c # Calculation frequency [Hz]

        self.tau = 1 / self.f_c  # seconds between state updates [s]

        self.cavity.simulation(self.k, self.f_c, self.E_in)
        
        self.f_calc = self.cavity.f_calc # Assigning to the self.f_c
        #print("Calculation frequency [Hz]: ", self.f_calc)

        self.t_a = self.cavity.t_a
        self.r_a = self.cavity.r_a
        self.r_b = self.cavity.r_b
        self.length = self.__L__
        self.F = self.cavity.Finesse()

        self.random_L = random_L
        self.displacement_factor = displacement_factor

        self.initial_length_deviation = 6. * self.lambd/(2.*self.cavity.Finesse())  # Multiplication factor is a magic number

        # STEP COUNTERS
        self.step_nmbr = 0

        # POWERS THRESHOLD
        self.p_thresh = 0.9999
        self.p_unlock_thresh = np.power(self.cavity.gain(), 2) * 0.5
        self.once_locked = False
        self.locked_steps = 0
        self.max_steps = max_steps

        # DECIDE IF THE ENVIRONMENT IS IN EVAL MODE
        self.eval_mode = False

        # NOISE
        self.P_noise_sigma = P_noise_sigma
        self.PDH_noise_sigma = PDH_noise_sigma

        # VELOCITY
        self.v = self.v_cr(self.cavity, self.lambd) # [m/s]
        self.v_max = v_max # Velocity factor 
        self.f_ref = 10e3 # Reference frequency for the velocity scaling [Hz]
        self.v_eff = self.v_max * (self.f_calc / self.f_ref) # Effective velocity factor
    
        # Lenght at which to fail the episode
        # self.lengh_limit = self.cavity.L + 3*lambd


        # SETTING THE ACTION AND OBSERVATION SPACES
        # self.action_box_limit = self.v_max/self.f_calc*self.v  # Maximum mirror displacement per step [m]
        self.action_box_limit = self.v_eff/self.f_calc*self.v
        max_noise_sigma = 3*max(P_noise_sigma, PDH_noise_sigma)
        self.action_space = spaces.Box(-1, 1, shape=(1,), dtype=np.float32)
        self.observation_space = spaces.Box(-1, 1, shape=(2,), dtype=np.float32)

        # DEFINE LISTS FOR POWERS, ACTIONS, REWARDS
        self.powers = []
        self.pdh = []
        if not self.no_render:
            self.actions = []
            self.rewards = []
            self.rewards_factors = {}
            self.score = 0.
            self.scores = []

        # REWARD FUNCTION
        self.reward_fn = get_reward_function(reward_fn)

        # TERMINATION CONDITION
        self.term_fn = get_termination_condition(term_fn)

        self.render_mode = render_mode
        
        self.state = None

    def step(self, actor_action):

        #self.resize_action_space(action)

        # Checking if action is within the action space and if it was called the reset function
        assert self.action_space.contains(
            actor_action
        ), f"{actor_action!r} ({type(actor_action)}) invalid"
        assert self.state is not None, "Call reset before using step method."

        action = (actor_action+1) * self.action_box_limit - self.action_box_limit
        self.length = self.length + action[0]

        # Storing previous state for reward calculation
        prev_state = self.state

        # Compute power for the choosen action and store it in the env state
        P, PDH = self.state
        E_out, _ = self.cavity.sim_step(d_zeta=action[0], E_in_laser=self.E_in)
        P = np.abs(E_out)**2
        P_normalized_pure = P / ((self.cavity.gain() * np.abs(self.E_in))**2 + 1e-12)

        if self.P_noise_sigma != 0.:
            P_normalized = P_normalized_pure + np.clip(np.random.normal(0., self.P_noise_sigma), -1, 1)
        else:
            P_normalized = P_normalized_pure
        gamma = 0.
        PDH = - np.angle(np.exp(gamma*1.j)*np.conjugate(self.E_in)*E_out)/np.pi

        # Normalization of Pound Drever Hall signal
        PDH = float(np.clip(PDH, -1.0, 1.0))

        if self.PDH_noise_sigma != 0.:
            PDH += np.clip(np.random.normal(0., self.PDH_noise_sigma), -self.PDH_noise_sigma*3, self.PDH_noise_sigma*3)
        
        self.powers.append(P_normalized)
        self.pdh.append(PDH)


        # DERIVATIVE OF THE POWER
        if len(self.powers) > 1:
            de = self.powers[-1] - self.powers[-2]
        else:
            de = 0.

        # COMPOSITION OF THE STATE
        self.state = (P_normalized, PDH)
    
        
        '''MAYBE TO DELETE'''
        # v_rel = abs(action[0]) / self.action_box_limit  # Normalized relative velocity
        # self.p_thresh = 0.9 + (0.9999 - 0.9) * np.exp(-5 * v_rel) # Threshold for the power to consider the cavity locked
	
	    # CHECKING IF THE CAVITY IS LOCKED
        if P_normalized > self.p_thresh and not self.once_locked:
            self.once_locked = True
            self.locked_steps = 1
        elif self.once_locked and P_normalized > self.p_thresh:
                self.locked_steps += 1
        else:
            self.locked_steps = 0
            self.once_locked = False
        #info = {"locked_steps": self.locked_steps, "once_locked": self.once_locked}
        # if self.locked_steps == 3:
        #     print("Locked at step: ", self.step_nmbr)
        #     print("The cavity length is {number_of_wavelength} half wavelengths".format(number_of_wavelength=self.length/(self.lambd/2.)))
        
                # CALCULATION OF THE REWARD
        reward_data = self.reward_fn(self.state, actor_action, info={"de": de, "action_box_limit": self.action_box_limit, 
                                                               "f_c": self.f_c, "powers": self.powers, 
                                                               "p_thresh": self.p_thresh, "once_locked": self.once_locked, 
                                                               "locked_steps": self.locked_steps, "step_nmbr": self.step_nmbr, 
                                                               "v_eff": self.v_eff, 'prev_obs': prev_state})
        reward = reward_data["total"]
        # Update the step counter and fill the state and action lists
        self.step_nmbr += 1
         
        # DEFINITION OF THE TERMINATION CONDITION
        info_term = {"step_nmbr": self.step_nmbr, "locked_steps": self.locked_steps, 
                "once_locked": self.once_locked, "p_thresh": self.p_thresh, "max_steps": self.max_steps}
        
        terminated, truncated = self.term_fn(self.state, action=action, info=info_term)

        # FILLING THE LISTS FOR RENDERING
        if not self.no_render:
            # rewards factors
            for key, value in reward_data.items():
                if key == "total":
                    continue
                if key not in self.rewards_factors:
                    self.rewards_factors[key] = []
                self.rewards_factors[key].append(value)

            self.rewards.append(reward)
            self.score+=reward
            self.scores.append(self.score)
            self.actions.append(actor_action)

        if self.render_mode == "human":
            self.render()

        #truncation=False as the time limit is handled by the `TimeLimit` wrapper added during `make`
        return np.array(self.state, dtype=np.float32), reward, terminated, truncated, info_term

    def reset(
        self,
        *,
        seed: Optional[int] = None,
        options: Optional[dict] = None,
    ):
        super().reset(seed=seed)
        
        # Reset the step number, power and lists
        self.step_nmbr = 0
        self.terminated = False
        self.locked_steps = 0

        P = 0.
        PDH = 0.
        info = {}
        self.powers=[]
        self.pdh = []

        self.once_locked = False

        if not self.no_render:
            self.actions=[]
            self.rewards=[]
            self.score = 0.
            self.scores=[]
            self.rewards_factors = {}

        # RESET THE CAVITY TO A RANDOM LENGTH
        if self.random_L:
            epsilon = np.random.uniform(self.initial_length_deviation*0.5, self.initial_length_deviation*2)
            epsilon_sign = np.random.uniform(-1., 1.)
            L = self.__L__ + epsilon * np.sign(epsilon_sign)*self.displacement_factor
            self.length = L
            # Reset the simulation
            if self.cavity_to_load is not None:
                self.cavity = self.cavity_to_load
                self.length = L
                self.cavity.__L__ = L
                self.cavity.simulation(self.k, self.f_c, self.E_in)
                self.f_calc = self.cavity.f_calc
                self.v_eff = self.v_max * (self.f_calc / self.f_ref) # Effective velocity factor
                self.action_box_limit = self.v_eff/self.f_calc*self.v

            else:
                self.cavity = op.Cavity(t_a=self.t_a , r_a=self.r_a , r_b=self.r_b , L=L)
                self.cavity.simulation(self.k, self.f_c, self.E_in)
                self.f_calc = self.cavity.f_calc
                self.v_eff = self.v_max * (self.f_calc / self.f_ref) # Effective velocity factor
                self.action_box_limit = self.v_eff/self.f_calc*self.v
        else:
            #self.cavity.sim_reset()
            if self.cavity_to_load is not None:
                self.cavity = self.cavity_to_load
                self.cavity.simulation(self.k, self.f_c, self.E_in)
                self.f_calc = self.cavity.f_calc
                self.v_eff = self.v_max * (self.f_calc / self.f_ref) # Effective velocity factor
                self.action_box_limit = self.v_eff/self.f_calc*self.v
            else:
                self.cavity = op.Cavity(t_a=self.t_a , r_a=self.r_a , r_b=self.r_b , L=self.__L__)
                self.cavity.simulation(self.k, self.f_c, self.E_in)
                self.f_calc = self.cavity.f_calc
                self.v_eff = self.v_max * (self.f_calc / self.f_ref) # Effective velocity factor
                self.action_box_limit = self.v_eff/self.f_calc*self.v

        self.state = (P, PDH)

        if self.render_mode == "human":
            self.render()

        return np.array(self.state, dtype=np.float32), info
    

    def render(self):
        if self.no_render:
            print("No render possible.")
            return
        
        if self.render_mode == "human":
            
            def animate(i, x, data):
                #d_zeta = env.action_space.sample()
                x += 1
                if self.terminated:
                    ani.event_source.stop()

                y = self.powers[i]
                data.append((x, y))
                ax.relim()
                ax.autoscale_view()
                line.set_data(*zip(*data))

            data = deque([(0, 0.)], maxlen=100)
            x=0
            y=0
            fig, ax = plt.subplots()
            
            line, = plt.plot(*zip(*data), c='black')

            ani = FuncAnimation(fig, animate, fargs=(x, data), frames=100, interval=10)
            plt.show()

        else:  
            fig, axs = plt.subplots(2, 2)
            plt.subplots_adjust(wspace=8.)
            axs[0, 0].grid(visible=True)
            axs[0, 0].plot(self.powers, label="cavity power")
            axs[0, 0].set_ylabel("Norm. cavity power")
            axs[0, 0].axhline(y=self.p_thresh, c='green', ls='--', lw='1')
            axs[0, 0].set_title("Input signals")

            ax2 = axs[0, 0].twinx()
            ax2.plot([0], label="cavity power")
            ax2.plot(self.pdh, label="PDH error")
            ax2.set_ylabel("PDH error")
            ax2.legend(loc=4)

            axs[1, 0].grid(visible=True)
            axs[1, 0].plot(self.actions)
            axs[1, 0].axhline(y=-self.action_box_limit, c='red', ls='--', lw='1')
            axs[1, 0].axhline(y=self.action_box_limit, c='red', ls='--', lw='1')
            axs[1, 0].set_title("Actions")
            axs[1, 0].sharex(axs[0, 0])
            axs[1, 0].set_ylabel("mirror displacement")
            axs[1, 0].set_xlabel("time step")

            axs[0, 1].grid(visible=True)
            axs[0, 1].plot(self.rewards, label='reward')
            axs[0, 1].set_ylabel("reward")
            axs[0, 1].set_title("Rewards")
            axs[0, 1].sharex(axs[0, 0])

            ax2 = axs[0, 1].twinx()
            ax2.plot([0], label="reward")
            for fact in self.rewards_factors:
                ax2.plot(self.rewards_factors[fact], label=fact, lw='0.7')
            ax2.set_ylabel("reward factors")
            ax2.legend()

            axs[1, 1].grid(visible=True)
            axs[1, 1].plot(self.scores)
            axs[1, 1].set_title("Scores")
            axs[1, 1].sharex(axs[0, 0])
            axs[1, 1].set_xlabel("time step")

            fig.tight_layout()
        
            return fig, axs
    
        return
    
    def close(self):
        plt.clf()


    # TO-DO
    # - Normalization of the action space between [0,1] (or [-1,1]?)
    # - Find a way to calculate self.steps_limit and the threshold (self.p_thresh) for the cavity power according to the velocity.
    # - Fix render function for live plot


    
