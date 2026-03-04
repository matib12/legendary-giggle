'''Simple integration of the limited exponential sinusoidal decay example using euler's method'''

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def f(t, y):
    return -2 * y * np.exp(-t) * np.sin(t)

# Time parameters
t0 = 0
tf = 5
dt = 0.1

# Time array
t = np.arange(t0, tf + dt, dt)

# Initialize solution array
y = np.zeros(len(t))

# Initial condition
y[0] = 1

# Euler's method
for i in range(1, len(t)):
    y[i] = y[i-1] + f(t[i-1], y[i-1]) * dt

# Plotting the results
# Compute exact solution using the integrated function
y_exact = np.exp(np.exp(-t)*(np.sin(t)+np.cos(t)) - 1)

# Create subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# First subplot: Euler method result
ax1.plot(t, y, label='Euler Method', marker='o')
ax1.set_title('Euler Method Integration')
ax1.set_xlabel('Time (t)')
ax1.set_ylabel('Solution (y)')
ax1.grid()
ax1.legend()

# Second subplot: Exact solution
ax2.plot(t, y_exact, label='Exact Solution', marker='s', color='orange')
ax2.set_title('Exact Solution')
ax2.set_xlabel('Time (t)')
ax2.set_ylabel('Solution (y)')
ax2.grid()
ax2.legend()

plt.tight_layout()
plt.savefig("output/simple_integration.png")
plt.show()
plt.close()
