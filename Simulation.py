import numpy as np
from Particle import Particle
import copy
import matplotlib.pyplot as plt

def print_particles(particles):
    for particle in particles:
        print(particle)

def plot_particle_properties(times, particle_properties, particles, property_name, ylabel):
    for i in range(len(particles)):
        plt.plot(times, [property_values[i] for property_values in particle_properties], label=f"{particles[i].name} - {property_name}")

    plt.xlabel('Time')
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()

def momentum_check(particles):
    total_momentum = 0
    for particle in particles:
        total_momentum += particle.momentum
    return total_momentum

def kinetic_energy(particles):
    total_kinetic_energy = 0
    for particle in particles:
        velocity_magnitude = np.linalg.norm(particle.velocity)
        kinetic_energy = 0.5 * particle.mass * velocity_magnitude**2
        total_kinetic_energy += kinetic_energy
    return total_kinetic_energy

# Read particle data from the input file
with open("input.txt", "r") as file:
    particle_data = file.read().strip().split('\n\n')  # Split particles by empty lines

# creates the particle in its class
particles = [Particle(data.split('\n')) for data in particle_data]

# Simulation using Runge-Kutta method
initial_gpe_rk4 = Particle.gpe(particles, particles)

Data_rk4 = []
time_rk4 = 0
initial_kinetic_energy_rk4 = kinetic_energy(particles)
initial_momentum_rk4 = 0
#runs the Runge-Kutta simulation for a set amount of time
for step in range(20000):
    for particle in particles:
        other_particles = [other_particle for other_particle in particles if other_particle != particle]
        if step == 0:
            initial_momentum_rk4 += particle.momentum
        particle.updateGravitationalAcceleration(other_particles)
    for particle in particles:
        particle.rk4update(1.0)
    time_rk4 += 1.0

    if step % 100 == 0:
        Data_rk4.append([time_rk4] + [copy.deepcopy(particle) for particle in particles])

# Save the generated data for Runge-Kutta
np.save("rk4method", Data_rk4, allow_pickle=True)

# Calculate final values for Runge-Kutta
final_momentum_rk4 = momentum_check(particles)
final_kinetic_energy_rk4 = kinetic_energy(particles)

# Print results for Runge-Kutta
print("Total Initial Momentum (RK4):", initial_momentum_rk4)
print("Total Momentum at the end (RK4):", final_momentum_rk4)
print("Percentage Difference in Momentum (RK4):", 100 * (initial_momentum_rk4 - final_momentum_rk4) / initial_momentum_rk4)

print("Initial Kinetic Energy (RK4):", initial_kinetic_energy_rk4)
print("Final Kinetic Energy (RK4):", final_kinetic_energy_rk4)
print("Percentage Difference in Kinetic Energy (RK4):", 100 * (initial_kinetic_energy_rk4 - final_kinetic_energy_rk4) / initial_kinetic_energy_rk4)

# Print gpe and final kinetic energy for RK4
final_gpe_rk4 = Particle.gpe(particles, particles)
initial_total_energy_rk4 = initial_kinetic_energy_rk4 + initial_gpe_rk4
final_total_energy_rk4 = final_kinetic_energy_rk4 + final_gpe_rk4

print("Initial Total Energy (RK4):", initial_total_energy_rk4)
print("Final Total Energy (RK4):", final_total_energy_rk4)

# Simulation using Euler method
initial_gpe_euler = Particle.gpe(particles, particles)

Data_euler = []
time_euler = 0
initial_kinetic_energy_euler = kinetic_energy(particles)
initial_momentum_euler = 0
#runs the eulers simulation for a set amount of time
for step in range(20000):
    for particle in particles:
        other_particles = [other_particle for other_particle in particles if other_particle != particle]
        if step == 0:
            initial_momentum_euler += particle.momentum
        particle.updateGravitationalAcceleration(other_particles)
    for particle in particles:
        particle.eulerupdate(1.0)
    time_euler += 1.0

    if step % 100 == 0:
        Data_euler.append([time_euler] + [copy.deepcopy(particle) for particle in particles])

# Save the generated data for Euler
np.save("Eulermethod", Data_euler, allow_pickle=True)

# Calculate final values for Euler and check momentum and energy is conserved
final_momentum_euler = momentum_check(particles)
final_kinetic_energy_euler = kinetic_energy(particles)

# Print results for Euler
print("Total Initial Momentum (Euler):", initial_momentum_euler)
print("Total Momentum at the end (Euler):", final_momentum_euler)
print("Percentage Difference in Momentum (Euler):", 100 * (initial_momentum_euler - final_momentum_euler) / initial_momentum_euler)

print("Initial Kinetic Energy (Euler):", initial_kinetic_energy_euler)
print("Final Kinetic Energy (Euler):", final_kinetic_energy_euler)
print("Percentage Difference in Kinetic Energy (Euler):", 100 * (initial_kinetic_energy_euler - final_kinetic_energy_euler) / initial_kinetic_energy_euler)

# Print gravitational potential energy and final kinetic energy for Euler
final_gpe_euler = Particle.gpe(particles, particles)
initial_total_energy_euler = initial_kinetic_energy_euler + initial_gpe_euler
final_total_energy_euler = final_kinetic_energy_euler + final_gpe_euler

print("Initial Total Energy (Euler):", initial_total_energy_euler)
print("Final Total Energy (Euler):", final_total_energy_euler)

# Plot graphs of properties for Runge-Kutta method
plot_particle_properties(
    [data[0] for data in Data_rk4],
    [[np.linalg.norm(getattr(particle, 'position')) for particle in data[1:]] for data in Data_rk4],
    particles,
    'Position',
    'Position for Runge-Kutta'
)

plot_particle_properties(
    [data[0] for data in Data_rk4],
    [[np.linalg.norm(getattr(particle, 'velocity')) for particle in data[1:]] for data in Data_rk4],
    particles,
    'Velocity',
    'Velocity for Runge-Kutta'
)

plot_particle_properties(
    [data[0] for data in Data_rk4],
    [[np.linalg.norm(getattr(particle, 'acceleration')) for particle in data[1:]] for data in Data_rk4],
    particles,
    'Acceleration',
    'Acceleration for Runge-Kutta'
)

# Plot graphs of properties for Euler
plot_particle_properties(
    [data[0] for data in Data_euler],
    [[np.linalg.norm(getattr(particle, 'position')) for particle in data[1:]] for data in Data_euler],
    particles,
    'Position',
    'Position for Euler'
)

plot_particle_properties(
    [data[0] for data in Data_euler],
    [[np.linalg.norm(getattr(particle, 'velocity')) for particle in data[1:]] for data in Data_euler],
    particles,
    'Velocity',
    'Velocity for Euler'
)

plot_particle_properties(
    [data[0] for data in Data_euler],
    [[np.linalg.norm(getattr(particle, 'acceleration')) for particle in data[1:]] for data in Data_euler],
    particles,
    'Acceleration',
    'Acceleration for Euler'
)
# When the graphs are generated, there is an option to save them; this is how I will import them to the report.
