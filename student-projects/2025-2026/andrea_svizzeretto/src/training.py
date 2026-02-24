'''
TRAINING SCRIPT

This script is part of TeRABIT Project.
The purpose of this script is to execute heavy training sessions for different RL agents on different simulated cavities.

Through some reward library functions it is possible to easily change the reward function used for training.

Partial observability and domain randomization are supported.

Some logging functionalities are included to keep track of the training progress.

Author: Andrea Svizzeretto
Date: January 2026
'''
# Import necessary libraries
import os
import datetime
from pyexpat import model
import time
import numpy as np
import gymnasium as gym
from gym_examples.envs.cavity import CavityEnv
import oreonspy as op
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from stable_baselines3 import DDPG, SAC, TD3
from sb3_contrib import RecurrentPPO
from stable_baselines3.common.noise import NormalActionNoise
from stable_baselines3.common.noise import OrnsteinUhlenbeckActionNoise
from stable_baselines3.common.callbacks import EvalCallback, StopTrainingOnRewardThreshold
from stable_baselines3.common.monitor import Monitor

def _format_table(data: dict):
    """
    Create a simple ASCII table with two columns:
    | Field | Value |
    """

    headers = ["Field", "Value"]
    rows = [(str(k), str(v)) for k, v in data.items()]

    # Compute column widths
    col_w = [
        max(len(headers[0]), max(len(r[0]) for r in rows)),
        max(len(headers[1]), max(len(r[1]) for r in rows)),
    ]

    def fmt_row(a, b):
        return f"| {a.ljust(col_w[0])} | {b.ljust(col_w[1])} |"

    sep = f"+-{'-' * col_w[0]}-+-{'-' * col_w[1]}-+"

    out = [sep, fmt_row(*headers), sep]
    for k, v in rows:
        out.append(fmt_row(k, v))
    out.append(sep)

    return "\n".join(out)

def extract_cavity_info(env):
    """
    Safely extract cavity-related parameters from the environment.
    """
    unwrapped = env.unwrapped

    # Cavity length
    L = getattr(unwrapped, "length", None)

    # v_max (agent velocity limit)
    v_max = getattr(unwrapped, "v_max", None)

    # Finesse
    finesse = getattr(unwrapped, "F", None)

    # r_a, r_b, t_a (mirror reflectivities)
    r_a = getattr(unwrapped, "r_a", None)
    r_b = getattr(unwrapped, "r_b", None)
    t_a = getattr(unwrapped, "t_a", None)


    return L, finesse, v_max, r_a, r_b, t_a

def write_run_report_txt(
    run_dir,
    model_name,
    timestamp,
    agent_type,
    reward_name,
    f_c,
    cavity_UUID,
    L,
    r_a,
    r_b,
    t_a,
    finesse,
    v_max,
    test_passed_bool,
    episode_rewards,
    lock_success_count,
    training_time,
    extra_info=None,
):
    """Write a txt report with a clear table summarizing the run."""

    episode_rewards = np.asarray(episode_rewards, dtype=float)
    mean_rew = float(np.mean(episode_rewards)) if len(episode_rewards) else float("nan")
    std_rew = float(np.std(episode_rewards)) if len(episode_rewards) else float("nan")

    headers = [
        "model_name",
        "agent_type",
        "reward_name",
        "f_c [Hz]",
        "cavity_UUID",
        "L [m]",
        "r_a",
        "r_b",
        "t_a",
        "Finesse",
        "v_max",
        "test_passed",
        "locks/episodes",
        "mean_reward(5ep)",
        "std_reward(5ep)",
        "training_time",
    ]

    rows = [[
        model_name,
        agent_type,
        reward_name,
        str(f_c),
        cavity_UUID,
        f"{L:.3f}" if L is not None else "N/A",
        f"{r_a:.5f}" if r_a is not None else "N/A",
        f"{r_b:.5f}" if r_b is not None else "N/A",
        f"{t_a:.5f}" if t_a is not None else "N/A",
        f"{finesse:.1f}" if finesse is not None else "N/A",
        f"{v_max:.3f}" if v_max is not None else "N/A",
        str(bool(test_passed_bool)),
        f"{lock_success_count}/5",
        f"{mean_rew:.6f}",
        f"{std_rew:.6f}",
        f"{training_time} s",
    ]]

    table = _format_table(dict(zip(headers, rows[0])))

    lines = []
    lines.append(f"Run timestamp: {timestamp}")
    lines.append(f"Run directory: {run_dir}")
    lines.append("")
    lines.append(table)
    lines.append("")

    if extra_info:
        lines.append("Extra info:")
        for k, v in extra_info.items():
            lines.append(f"- {k}: {v}")

    report_name = f"{model_name}_{timestamp}_test_passed:{bool(test_passed_bool)}.txt"
    report_path = os.path.join(run_dir, report_name)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return report_path

def train_and_test_agent(env, eval_env, agent_type, reward_name, total_timesteps, log_dir, f_c, cavity_UUID, episode = 5, seed=0):
    """Function to train an RL agent on a specified environment with a given reward function."""
    
    # Create the environment
    #env = gym.make(env_name)
    
    # Set random seed for reproducibility
    env.reset(seed=seed)
    env.action_space.seed(seed)
    np.random.seed(seed)

    callback_on_best = StopTrainingOnRewardThreshold(reward_threshold=3000, verbose=1)
    eval_callback = EvalCallback(env, callback_on_new_best=callback_on_best, verbose=1, eval_freq=20000)

    # Define action noise for exploration (for DDPG and TD3)
    n_actions = env.action_space.shape[-1]
    action_noise = NormalActionNoise(mean=np.zeros(n_actions), sigma=0.1 * np.ones(n_actions))
    action_noise = OrnsteinUhlenbeckActionNoise(
        mean=np.zeros(n_actions),
        sigma=0.3 * np.ones(n_actions),
        theta=0.15
    )
    # --- Build model name and initialize the agent ---
    # The "model_name" is used for folder/file names as requested.
    model_name = f"{agent_type}_{reward_name}_fc{f_c}_cavity{cavity_UUID}"

    # Initialize the agent
    if agent_type == 'DDPG':
        model = DDPG('MlpPolicy', env, action_noise=action_noise, verbose=0, tensorboard_log=log_dir)
    elif agent_type == 'TD3':
        model = TD3('MlpPolicy', env, action_noise=action_noise, verbose=0, tensorboard_log=log_dir)
    elif agent_type == 'SAC':
        model = SAC('MlpPolicy', env, verbose=0, tensorboard_log=log_dir)
    elif agent_type == 'RecurrentPPO':
        model = RecurrentPPO('MlpLstmPolicy', env, verbose=0, tensorboard_log=log_dir)
    else:
        raise ValueError("Unsupported agent type. Choose from 'DDPG', 'TD3', or 'SAC'.\n\n")
    
    # Train the agent
    print("\n\n----------------------TRAINING PHASE...----------------------------\n")
    start_time = time.time()
    model.learn(total_timesteps=total_timesteps, callback=eval_callback)
    end_time = time.time()
    training_dt = end_time - start_time
    print(f"Training for {agent_type}_{reward_name}_{f_c} completed. | Training time: {training_dt} |\n\n")
    model_path = os.path.join(log_dir, f"{agent_type}_{reward_name}_model")
    model.save(model_path)
    print(f"- Model saved at {model_path}. -\n\n\n")
    del model
    # Load the trained agent
    if agent_type == 'DDPG':
        model = DDPG.load(model_path, env=eval_env)
    elif agent_type == 'TD3':
        model = TD3.load(model_path, env=eval_env)
    elif agent_type == 'SAC':
        model = SAC.load(model_path, env=eval_env)
    elif agent_type == 'RecurrentPPO':
        model = RecurrentPPO.load(model_path, env=eval_env)
    else:
        raise ValueError("Model can't be loaded\n")
    # Test the trained agent
    print("----------------------TESTING PHASE----------------------------\n\n")
    states = []
    locked_steps = 0.
    
    episode_rewards = []
    lock_success_count = 0
    hold_success_steps = 20
    test_passed_bool = False

    for ep in range(1, episode+1):
        obs, _ = eval_env.reset()
        done = False
        score = 0
        lstm_states = None
        n_envs = getattr(eval_env, 'num_envs', 1)
        episode_start = np.ones((1,), dtype=bool)
        terminated = False
        if agent_type == 'RecurrentPPO':
            while not done:
                d_zeta, lstm_states = model.predict(obs, state = lstm_states, episode_start=episode_start, deterministic=True)
                obs, reward, terminated, truncated, info = eval_env.step(np.array(d_zeta, dtype=np.float32))
                states.append(obs)
                score+=reward
                # last_info = info  # Store the last info dictionary
                last_terminated = terminated
                done = bool(terminated or truncated)
                episode_start[:] = done
        else:
            while not done:
                d_zeta, _states = model.predict(obs, deterministic=True)
                obs, reward, terminated, truncated, info = eval_env.step(np.array(d_zeta, dtype=np.float32))
                states.append(obs)
                score+=reward
                last_terminated = terminated
                # last_info = info  # Store the last info dictionary
                done = bool(terminated or truncated)
        #locked_steps = int(env.unwrapped.locked_steps)
        episode_rewards.append(score)
        # Determine if lock happened:
        # locked_steps = int(last_info.get("locked_steps", 0))

        # Fallback threshold (in case termination_reason is missing)
        ep_locked = bool(locked_steps >= hold_success_steps)

        # if ep_locked:
        #     lock_success_count += 1
        if last_terminated:
            lock_success_count += 1
        print(f"-----------------------------EPISODE: {ep}/{episode}-------------------------------")
        print(f"- Test {ep}/{episode} on {model_name} -> Score: {score:.3f} | locked={ep_locked} -\n") 
        fig, axs = env.render()
        fig.savefig(os.path.join(log_dir, f"{agent_type}_{reward_name}_test_Episode_{ep}.pdf"))
        plt.close(fig)
        print(f"- Test plot saved at {log_dir}/{agent_type}_{reward_name}_test_Episode_{ep}.pdf -\n")
        print(f"- Test on {agent_type}_{reward_name}_{f_c}_{cavity_UUID} -> Score:{score} -\n\n")
        
        if lock_success_count >= 1:
            test_passed_bool = True
            print(f"- Test PASSED: locked episodes = {lock_success_count} /5 -\n")
        else:
            print(f"- Test FAILED: locked_episodes = {lock_success_count} /5 -\n")
    print("------------------------------------------------------------\n\n")
    return model_name, model_path, test_passed_bool, episode_rewards, lock_success_count, training_dt

base_dir = "./optical_cavities_testset/"
file_names = os.listdir(base_dir)
file_paths = [os.path.join(base_dir, file_name) for file_name in file_names]

# Extract UUIDs from file names in tests/optical_cavities_testset
selected_cavities = []
selected_cavities_UUIDS = []
filtered_cavities = []
filtered_UUIDS = []

for file_path in file_paths:
    file_name = os.path.basename(file_path)

    if not file_name.endswith(".xml"):
        continue

    name_no_ext = os.path.splitext(file_name)[0]
    uuid = name_no_ext.split("_")[1]

    cavity = op.Cavity()
    cavity.xml_load(file_path)

    selected_cavities_UUIDS.append(uuid)
    selected_cavities.append(cavity)

print(f"Selected {len(selected_cavities)} cavities.") 

# Select one cavity every five
filtered_cavities = selected_cavities[::5]
filtered_UUIDS = selected_cavities_UUIDS[::5]
print(f"Filtered {len(filtered_cavities)} cavities from the selected list.")


# -----------------------
# TRAINING LOOP
# -----------------------
rewards = ["NewReward_ActionAware", 
           "ModifiedReward_ActionAware", 
           "RealEnv", 
           "LockSteps_v1_Modified_NewReward_ActionAware", 
           "LockSteps_v2_Modified_NewReward_ActionAware", 
           "LockSteps_v3_Modified_NewReward_ActionAware",
           'reward_stable',
           'stepsize_RealEnv',
           'curiosity_lock_reward',
           'fine_curiosity_lock_reward',
           'fine_curiosity_lock_reward_v2']
agents = ["DDPG", "TD3", "SAC", "RecurrentPPO"]
f_c_list = [200, 1000, 5000, 10000, 150e3]
total_timesteps = 200000
max_steps = [200, 300, 500, 1000, 2000, 4000]

# Root log folder as requested
training_logs_root = os.path.join("training", "logs")
os.makedirs(training_logs_root, exist_ok=True)


for reward_name in rewards:
    for agent_type in agents:
        if reward_name == 'fine_curiosity_lock_reward_v2':
            print(f"\n\n================= TRAINING {agent_type} with reward: {reward_name} =================\n\n")
        for cavity, uuid in zip(filtered_cavities, filtered_UUIDS):
            for f_c in f_c_list:
                for max_step in max_steps:
                    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
                    # Single run directory:
                    # training/logs/{model_name}/{timestamp}/
                    # NOTE: model_name is constructed inside train_and_test_agent,
                    # but we need the folder before training to pass it to tensorboard_log.
                    # So we build it here with the same naming convention.
                    model_name_preview = f"{agent_type}_{reward_name}_fc-{f_c}Hz_cavity-{uuid}_maxstep{max_step}"
                    run_dir = os.path.join(training_logs_root, model_name_preview, f"{timestamp}_{agent_type}")
                    os.makedirs(run_dir, exist_ok=True)
            
            
                    env = gym.make(
                        "gym_examples/Cavity-v0",
                        cavity_to_load=cavity,
                        reward_fn=reward_name,  
                        t_a=0.01377,
                        r_a=0.986,
                        r_b=0.99999,
                        L=3000.0,
                        f_c=f_c,
                        v_max=5.0,
                        max_steps = max_step,
                        random_L=True,
                        in_resonance_L=False,
                        displacement_factor=1.0,
                    )
                    eval_env = gym.make(
                        "gym_examples/Cavity-v0",
                        cavity_to_load=cavity,
                        reward_fn=reward_name,  
                        t_a=0.01377,
                        r_a=0.986,
                        r_b=0.99999,
                        L=3000.0,
                        f_c=f_c,
                        v_max=5.0,
                        random_L=True,
                        max_steps = max_step,
                        in_resonance_L=True,
                        displacement_factor=3.0,
                    )
            
                    env = Monitor(env)  # Wrap the environment with Monitor for logging
                    eval_env = Monitor(eval_env)  # Wrap the evaluation environment as well

                    model_name, model_path, test_passed_bool, episode_rewards, lock_success_count, training_dt = train_and_test_agent(
                        env=env,
                        eval_env=eval_env,
                        agent_type=agent_type,
                        reward_name=reward_name,
                        total_timesteps=total_timesteps,
                        log_dir=run_dir,
                        f_c=f_c,
                        cavity_UUID=uuid,
                        episode=5,
                        seed=42,
                    )
                    # Extract cavity/environment parameters
                    L, finesse, v_max, r_a, r_b, t_a = extract_cavity_info(env)

                    report_path = write_run_report_txt(
                        run_dir=run_dir,
                        model_name=model_name,
                        timestamp=timestamp,
                        agent_type=agent_type,
                        reward_name=reward_name,
                        f_c=f_c,
                        cavity_UUID=uuid,
                        L=L,
                        r_a=r_a,
                        r_b=r_b,
                        t_a=t_a,
                        finesse=finesse,
                        v_max=v_max,
                        test_passed_bool=test_passed_bool,
                        episode_rewards=episode_rewards,
                        lock_success_count=lock_success_count,
                        training_time=training_dt,
                        extra_info={
                            "model_path": model_path,
                            "total_timesteps": total_timesteps,
                            "max_steps": max_step,
                        }
                    )
                    print(f"\n\nReport saved: {report_path}\n\n")
                    env.close()
                    eval_env.close()
