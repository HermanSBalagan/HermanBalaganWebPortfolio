import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

# Landing Gear Load Calculation (Aerodynamic + Ground Reaction Forces)
def landing_gear_load(aircraft_weight, impact_load_factor, num_wheels):
    """Calculates the force on each landing gear wheel during touchdown."""
    total_load = aircraft_weight * impact_load_factor
    wheel_load = total_load / num_wheels
    return wheel_load

# Hydraulic Actuator Force Calculation (Pascal’s Law)
def hydraulic_actuator_force(pressure, piston_diameter):
    """Calculates force exerted by hydraulic actuators in the landing gear."""
    area = np.pi * (piston_diameter / 2) ** 2
    force = pressure * area
    return force

# Shock Absorber Energy Dissipation Calculation
def shock_absorber_energy(mass, velocity, stroke_length):
    """Estimates energy absorbed by the shock absorber during landing."""
    kinetic_energy = 0.5 * mass * velocity ** 2
    absorber_efficiency = 0.85  # Assumed efficiency
    absorbed_energy = kinetic_energy * absorber_efficiency
    required_force = absorbed_energy / stroke_length
    return absorbed_energy, required_force

# Tire Pressure Calculation (Load vs. Tire Contact Patch Area)
def tire_pressure(wheel_load, contact_area):
    """Calculates the required tire pressure based on wheel load and contact patch area."""
    return wheel_load / contact_area

# Monte Carlo Simulation for Failure Probability of Landing Gear Components
def monte_carlo_failure(prob_failure, n_simulations=10000):
    """Simulates failure occurrences based on probability."""
    failures = np.random.rand(n_simulations) < prob_failure
    return np.sum(failures) / n_simulations

# Recommended Maintenance Interval Calculation
def maintenance_interval(failure_rate, flight_cycles):
    """Estimates recommended time between maintenance based on failure rate and usage."""
    if failure_rate == 0:
        return float('inf')  # No maintenance required if failure rate is zero
    return flight_cycles / failure_rate

# Example Calculation for Aircraft Landing Gear
if __name__ == "__main__":
    # Example Parameters
    aircraft_weight = 50000  # in Newtons
    impact_load_factor = 1.5  # Load multiplier during landing impact
    num_wheels = 4  # Total wheels in the landing gear system
    hydraulic_pressure = 20000000  # Pascal
    piston_diameter = 0.15  # meters
    mass = 5000  # kg per landing gear strut
    velocity = 3  # m/s (touchdown vertical velocity)
    stroke_length = 0.5  # meters
    contact_area = 0.2  # m² per tire
    prob_failure = 0.02  # 2% Failure Probability
    flight_cycles = 10000  # Number of flight cycles before expected failure
    
    # Compute Landing Gear Load
    wheel_load = landing_gear_load(aircraft_weight, impact_load_factor, num_wheels)
    print(f"Load per Wheel: {wheel_load:.2f} N")
    
    # Compute Hydraulic Actuator Force
    actuator_force = hydraulic_actuator_force(hydraulic_pressure, piston_diameter)
    print(f"Hydraulic Actuator Force: {actuator_force:.2f} N")
    
    # Compute Shock Absorber Energy Dissipation
    absorbed_energy, required_force = shock_absorber_energy(mass, velocity, stroke_length)
    print(f"Shock Absorber Energy Dissipation: {absorbed_energy:.2f} J")
    print(f"Required Force on Shock Absorber: {required_force:.2f} N")
    
    # Compute Required Tire Pressure
    pressure = tire_pressure(wheel_load, contact_area)
    print(f"Required Tire Pressure: {pressure:.2f} Pa")
    
    # Monte Carlo Failure Simulation
    failure_rate = monte_carlo_failure(prob_failure)
    print(f"Predicted Failure Rate: {failure_rate * 100:.2f}%")
    
    # Compute Recommended Maintenance Interval
    recommended_maintenance = maintenance_interval(failure_rate, flight_cycles)
    print(f"Recommended Maintenance Interval: {recommended_maintenance:.2f} flight cycles")
    
    # Visualization of Failure Probability
    probabilities = np.linspace(0, 0.1, 50)
    failure_rates = [monte_carlo_failure(p) for p in probabilities]
    
    plt.plot(probabilities * 100, np.array(failure_rates) * 100, label="Failure Probability")
    plt.xlabel("Component Failure Probability (%)")
    plt.ylabel("Simulated Failure Rate (%)")
    plt.title("Monte Carlo Failure Prediction for Landing Gear")
    plt.legend()
    plt.grid()
    plt.show()
