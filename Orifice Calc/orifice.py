import math
import CoolProp.CoolProp as CP

def calculate_compressible_orifice_flow_rate(C_d, A, P1, P2, fluid, T1, gamma=None):
    """
    Calculate the mass flow rate through a compressible orifice using CoolProp for fluid properties.

    Parameters:
    - C_d: Discharge coefficient (dimensionless)
    - A: Orifice area (m²)
    - P1: Inlet pressure (Pa)
    - P2: Outlet pressure (Pa)
    - fluid: Fluid name (string, e.g., "Air")
    - T1: Inlet temperature (K)
    - gamma: Ratio of specific heats (optional, CoolProp can provide it)

    Returns:
    - Mass flow rate (dot_m) in kg/s
    """
    # Use CoolProp to get the density at inlet
    rho = CP.PropsSI('D', 'P', P1, 'T', T1, fluid)  # Density at inlet (kg/m³)

    # Calculate gamma (specific heat ratio)
    if gamma is None:
        Cp = CP.PropsSI('C', 'P', P1, 'T', T1, fluid)  # Specific heat at constant pressure (J/kg·K)
        Cv = CP.PropsSI('O', 'P', P1, 'T', T1, fluid)  # Specific heat at constant volume (J/kg·K)
        R = Cp - Cv  # Specific gas constant (J/kg·K)
        gamma = Cp / Cv  # Calculate specific heat ratio

    # Check if the flow is choked
    if P2 / P1 < (2 / (gamma + 1))**(gamma / (gamma - 1)):  # Choked flow condition
        mdot = C_d * A * P1 * math.sqrt(gamma / R / T1 * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1)))
        return mdot
    else:
        # Subsonic flow
        mdot = C_d * A * rho * (P2 / P1) ** (1 / gamma) * math.sqrt(2 * R * T1 * (gamma / (gamma - 1)) * (1 - (P2 / P1) ** ((gamma - 1) / gamma)))
        return mdot


def calculate_orifice_exit_velocity(C_d, P1, P2, rho, gamma):
    """
    Calculate the exit velocity of the gas through a compressible orifice.

    Parameters:
    - C_d: Discharge coefficient (dimensionless)
    - P1: Inlet pressure (Pa)
    - P2: Outlet pressure (Pa)
    - rho: Fluid density at inlet (kg/m³)
    - gamma: Ratio of specific heats (dimensionless)

    Returns:
    - Exit velocity (v) in m/s
    """
    if P2 / P1 < (2 / (gamma + 1))**(gamma / (gamma - 1)):  # Choked flow condition
        # Choked flow velocity (sonic velocity)
        return math.sqrt(gamma * (2 / (gamma + 1)) * (P1 / rho) * ((gamma - 1) / (gamma + 1)))
    else:
        # Subsonic flow velocity
        return math.sqrt(2 * (P1 - P2) / rho)


def main():
    # Input parameters
    C_d = 0.61  # Discharge coefficient (example value for sharp-edged orifice)
    diameter = 0.05  # Diameter of the orifice (m)
    P1 = 500000  # Inlet pressure (Pa)
    P2 = 100000  # Outlet pressure (Pa)
    fluid = "Helium"  # Fluid name (e.g., "Air", "Water", etc.)
    T1 = 300  # Inlet temperature (K)
    
    # Calculate the orifice area (A = π * (d/2)^2)
    A = math.pi * (diameter / 2)**2

    # Calculate the mass flow rate
    mass_flow_rate = calculate_compressible_orifice_flow_rate(C_d, A, P1, P2, fluid, T1)
    print(f"Mass flow rate through the orifice: {mass_flow_rate:.4f} kg/s")

    # Calculate the density and specific heat ratio using CoolProp
    rho = CP.PropsSI('D', 'P', P1, 'T', T1, fluid)  # Density (kg/m³)
    print(f"Upstream density: {rho:.4f} kg/m^3")

    # Calculate the specific heat ratio gamma
    Cp = CP.PropsSI('C', 'P', P1, 'T', T1, fluid)  # Specific heat at constant pressure (J/kg·K)
    print(f"Upstream C_p: {Cp:.4f} J/kg/K")
    Cv = CP.PropsSI('O', 'P', P1, 'T', T1, fluid)  # Specific heat at constant volume (J/kg·K)
    print(f"Upstream C_v: {Cv:.4f} J/kg/K")
    R = Cp - Cv  # Specific gas constant (J/kg·K)
    print(f"Upstream R: {R:.4f} J/kg/K")
    gamma = Cp / Cv  # Specific heat ratio
    print(f"Upstream gamma: {gamma:.4f}")

    # Calculate the exit velocity of the gas
    velocity = calculate_orifice_exit_velocity(C_d, P1, P2, rho, gamma)
    print(f"Exit velocity of the gas: {velocity:.4f} m/s")


if __name__ == "__main__":
    main()
