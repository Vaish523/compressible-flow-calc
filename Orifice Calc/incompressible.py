import math
import CoolProp.CoolProp as CP

def calculate_incompressible_orifice_flow_rate(C_d, A, P1, P2, fluid, T1):
    """
    Calculate the mass flow rate through an incompressible orifice using CoolProp for fluid properties.

    Parameters:
    - C_d: Discharge coefficient (dimensionless)
    - A: Orifice area (m²)
    - P1: Inlet pressure (Pa)
    - P2: Outlet pressure (Pa)
    - fluid: Fluid name (string, e.g., "Water")
    - T1: Inlet temperature (K)

    Returns:
    - Mass flow rate (dot_m) in kg/s
    - Exit velocity (v) in m/s
    """
    # Calculate the pressure drop
    delta_P = P1 - P2

    # Use CoolProp to get the density at the inlet (in kg/m³)
    rho = CP.PropsSI('D', 'P', P1, 'T', T1, fluid)

    # Calculate the mass flow rate (in kg/s) using the incompressible orifice formula
    mass_flow_rate = C_d * A * math.sqrt(2 * delta_P * rho)

    # Calculate the exit velocity (in m/s)
    velocity = math.sqrt(2 * delta_P / rho)

    return mass_flow_rate, velocity

def main():
    # Input parameters
    C_d = 0.61  # Discharge coefficient (dimensionless)
    diameter = 0.05  # Diameter of the orifice (m)
    P1 = 500000  # Inlet pressure (Pa)
    P2 = 100000  # Outlet pressure (Pa)
    fluid = "Water"  # Fluid name (e.g., "Water", "Ethylene glycol", etc.)
    T1 = 300  # Inlet temperature (K)

    # Calculate the orifice area (A = π * (d/2)^2)
    A = math.pi * (diameter / 2)**2

    # Calculate the mass flow rate and exit velocity
    mass_flow_rate, velocity = calculate_incompressible_orifice_flow_rate(C_d, A, P1, P2, fluid, T1)

    # Output the results
    print(f"Mass flow rate through the orifice: {mass_flow_rate:.4f} kg/s")
    print(f"Exit velocity of the fluid: {velocity:.4f} m/s")

if __name__ == "__main__":
    main()
