import numpy as np
import matplotlib.pyplot as plt

length = np.pi
num_points = 10
dx = length / num_points
dt = 0.2
time_steps = 50

displacement = np.array([np.sin(i * dx) for i in range(num_points + 1)])
velocity = np.zeros(num_points + 1)

kinetic_energy_list = []
potential_energy_list = []
total_energy_list = []

def compute_acceleration(u, dx):
    acceleration = np.zeros_like(u)
    for i in range(1, num_points):
        acceleration[i] = (u[i - 1] - 2 * u[i] + u[i + 1]) / dx**2
    return acceleration

# midpoint calculations
def simulate_step(u, v, a, dt, dx):
    half_step_velocity = v + 0.5 * a * dt
    new_displacement = u + half_step_velocity * dt
    new_displacement[0] = 0
    new_displacement[-1] = 0
    new_acceleration = compute_acceleration(new_displacement, dx)
    new_velocity = half_step_velocity + 0.5 * new_acceleration * dt
    return new_displacement, new_velocity, new_acceleration

def compute_energy(u, v, dx):
    kinetic = 0.0
    potential = 0.0
    for i in range(1, num_points):
        kinetic += 0.5 * dx * v[i]**2
        potential += 0.5 * ((u[i + 1] - u[i])**2) / dx
    return kinetic, potential

acceleration = compute_acceleration(displacement, dx)

for _ in range(time_steps):
    displacement, velocity, acceleration = simulate_step(displacement, velocity, acceleration, dt, dx)
    kin_energy, pot_energy = compute_energy(displacement, velocity, dx)
    total = kin_energy + pot_energy
    kinetic_energy_list.append(kin_energy)
    potential_energy_list.append(pot_energy)
    total_energy_list.append(total)

time_array = np.arange(time_steps) * dt

plt.figure(figsize=(10, 6))
plt.plot(time_array, kinetic_energy_list, label='Kinetic Energy')
plt.plot(time_array, potential_energy_list, label='Potential Energy')
plt.plot(time_array, total_energy_list, label='Total Energy')
plt.xlabel('Time [s]')
plt.ylabel('Energy')
plt.title('Energy Evolution of Vibrating String')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
