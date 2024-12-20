import math
from CoolProp.CoolProp import PropsSI

def calculate_venturi_parameters(chamber_pressure, dp_injector, line_losses, cavitation_ratio, Cd, thrust, isp, of_ratio, fluid, temperature, tube_od_fuel, tube_wall_fuel):
    """
    Calculate venturi parameters for a liquid rocket engine Fuel line.

    Parameters:
    - chamber_pressure: Expected chamber pressure (psi)
    - dp_injector: Pressure drop across injector to chamber (psi)
    - line_losses: Assumed line losses between venturi exit and injector inlet (psi)
    - cavitation_ratio: Venturi pressure ratio needed for cavitation (e.g., 0.85)
    - Cd: Discharge coefficient of the venturi (dimensionless)
    - thrust: Required thrust (N)
    - isp: Specific impulse (s)
    - of_ratio: Oxidizer-to-fuel mass ratio (dimensionless)
    - fluid: Name of the fluid (e.g., 'nDodecane' for RP-1)
    - temperature: Fluid temperature (K)
    - tube_od_fuel: Outer diameter of the fuel tube (inches)
    - tube_wall_fuel: Wall thickness of the fuel tube (inches)

    Returns:
    - venturi_exit_pressure: Pressure at venturi exit (psi)
    - venturi_inlet_pressure: Pressure at venturi inlet (psi)
    - required_CdA: Required CdA (m^2)
    - venturi_throat_diameter: Throat diameter (inches)
    - tube_inner_diameter_fuel: Inner diameter of the fuel tube (inches)
    """
    # Calculate total mass flow rate from thrust and Isp
    g0 = 9.80665  # Standard gravity (m/s^2)
    total_flow_rate = thrust / (isp * g0)  # Total mass flow rate (kg/s)

    # Calculate oxidizer and fuel flow rates
    oxidizer_flow_rate = total_flow_rate * (of_ratio / (1 + of_ratio))
    fuel_flow_rate = total_flow_rate - oxidizer_flow_rate

    # Convert pressures from psi to Pa
    chamber_pressure_pa = chamber_pressure * 6894.76
    dp_injector_pa = dp_injector * 6894.76
    line_losses_pa = line_losses * 6894.76

    # Calculate venturi exit pressure in Pa
    venturi_exit_pressure_pa = chamber_pressure_pa + dp_injector_pa + line_losses_pa

    # Calculate venturi inlet pressure using cavitation pressure ratio
    venturi_inlet_pressure_pa = venturi_exit_pressure_pa / cavitation_ratio

    # Get fluid density and vapor pressure at inlet conditions
    density = PropsSI('D', 'P', venturi_inlet_pressure_pa, 'T', temperature, fluid)  # kg/m^3
    vapor_pressure = PropsSI('P', 'Q', 0, 'T', temperature, fluid)  # Pa
    print(vapor_pressure)

    # Required CdA (discharge coefficient * area)
    required_CdA = fuel_flow_rate / (Cd * math.sqrt(2 * density * (venturi_inlet_pressure_pa - vapor_pressure)))

    # Throat area
    throat_area = required_CdA / Cd

    # Throat diameter (from area)
    venturi_throat_diameter_m = math.sqrt((4 * throat_area) / math.pi)

    # Convert pressures back to psi
    venturi_exit_pressure = venturi_exit_pressure_pa / 6894.76
    venturi_inlet_pressure = venturi_inlet_pressure_pa / 6894.76

    # Calculate inner diameter of the fuel tube
    tube_od_fuel_m = tube_od_fuel * 0.0254  # Convert inches to meters
    tube_wall_fuel_m = tube_wall_fuel * 0.0254  # Convert inches to meters
    tube_inner_diameter_fuel_m = tube_od_fuel_m - 2 * tube_wall_fuel_m

    # Convert diameters to inches
    venturi_throat_diameter = venturi_throat_diameter_m / 0.0254
    tube_inner_diameter_fuel = tube_inner_diameter_fuel_m / 0.0254

    return venturi_exit_pressure, venturi_inlet_pressure, required_CdA, venturi_throat_diameter, tube_inner_diameter_fuel

def calculate_oxidizer_venturi(chamber_pressure, dp_injector, line_losses, cavitation_ratio, Cd, oxidizer_flow_rate, fluid, temperature, tube_od_oxidizer, tube_wall_oxidizer):
    """
    Calculate venturi parameters for a liquid oxygen (LOX) line.

    Parameters:
    - chamber_pressure: Expected chamber pressure (psi)
    - dp_injector: Pressure drop across injector to chamber (psi)
    - line_losses: Assumed line losses between venturi exit and injector inlet (psi)
    - cavitation_ratio: Venturi pressure ratio needed for cavitation (e.g., 0.85)
    - Cd: Discharge coefficient of the venturi (dimensionless)
    - oxidizer_flow_rate: Mass flow rate of oxidizer (kg/s)
    - fluid: Name of the fluid (e.g., 'LOX' for liquid oxygen)
    - temperature: Fluid temperature (K)
    - tube_od_oxidizer: Outer diameter of the oxidizer tube (inches)
    - tube_wall_oxidizer: Wall thickness of the oxidizer tube (inches)

    Returns:
    - venturi_exit_pressure: Pressure at venturi exit (psi)
    - venturi_inlet_pressure: Pressure at venturi inlet (psi)
    - required_CdA: Required CdA (m^2)
    - venturi_throat_diameter: Throat diameter (inches)
    - tube_inner_diameter_oxidizer: Inner diameter of the oxidizer tube (inches)
    """
    # Convert pressures from psi to Pa
    chamber_pressure_pa = chamber_pressure * 6894.76
    dp_injector_pa = dp_injector * 6894.76
    line_losses_pa = line_losses * 6894.76

    # Calculate venturi exit pressure in Pa
    venturi_exit_pressure_pa = chamber_pressure_pa + dp_injector_pa + line_losses_pa

    # Calculate venturi inlet pressure using cavitation pressure ratio
    venturi_inlet_pressure_pa = venturi_exit_pressure_pa / cavitation_ratio

    # Get fluid density and vapor pressure at inlet conditions
    density = PropsSI('D', 'P', venturi_inlet_pressure_pa, 'T', temperature, fluid)  # kg/m^3
    vapor_pressure = PropsSI('P', 'Q', 0, 'T', temperature, fluid)  # Pa
    print(vapor_pressure)

    # Required CdA (discharge coefficient * area)
    required_CdA = oxidizer_flow_rate / (Cd * math.sqrt(2 * density * (venturi_inlet_pressure_pa - vapor_pressure)))

    # Throat area
    throat_area = required_CdA / Cd

    # Throat diameter (from area)
    venturi_throat_diameter_m = math.sqrt((4 * throat_area) / math.pi)

    # Convert pressures back to psi
    venturi_exit_pressure = venturi_exit_pressure_pa / 6894.76
    venturi_inlet_pressure = venturi_inlet_pressure_pa / 6894.76

    # Calculate inner diameter of the oxidizer tube
    tube_od_oxidizer_m = tube_od_oxidizer * 0.0254  # Convert inches to meters
    tube_wall_oxidizer_m = tube_wall_oxidizer * 0.0254  # Convert inches to meters
    tube_inner_diameter_oxidizer_m = tube_od_oxidizer_m - 2 * tube_wall_oxidizer_m

    # Convert diameters to inches
    venturi_throat_diameter = venturi_throat_diameter_m / 0.0254
    tube_inner_diameter_oxidizer = tube_inner_diameter_oxidizer_m / 0.0254

    return venturi_exit_pressure, venturi_inlet_pressure, required_CdA, venturi_throat_diameter, tube_inner_diameter_oxidizer


# Example Inputs
chamber_pressure = 500.0  # Chamber pressure in psi (e.g., 70 bar)
dp_injector = 0.2 * chamber_pressure  # Injector pressure drop (20% of chamber pressure)
line_losses = 0.05 * chamber_pressure  # Line losses (5% of chamber pressure)
cavitation_ratio = 0.85  # Venturi pressure ratio for cavitation
Cd = 0.95  # Discharge coefficient
thrust = 4448.0  # Required thrust in N
isp = 346  # Specific impulse in seconds
of_ratio = 3.9  # Oxidizer-to-fuel ratio
fuel_temperature = 110  # temperature in K
oxidizer_temperature = 90  # temperature in K
fuel_fluid = 'Methane'  # RP-1 surrogate fluid
oxidizer_fluid = 'Oxygen'  # Liquid oxygen
tube_od_fuel = 0.5  # Fuel tube outer diameter in inches
tube_wall_fuel = 0.049  # Fuel tube wall thickness in inches
tube_od_oxidizer = 0.5  # Oxidizer tube outer diameter in inches
tube_wall_oxidizer = 0.049  # Oxidizer tube wall thickness in inches

# Calculate RP-1 venturi parameters
venturi_exit_pressure, venturi_inlet_pressure, required_CdA, venturi_throat_diameter, tube_inner_diameter_fuel = calculate_venturi_parameters(
    chamber_pressure, dp_injector, line_losses, cavitation_ratio, Cd, thrust, isp, of_ratio, fuel_fluid, fuel_temperature, tube_od_fuel, tube_wall_fuel
)

# Calculate oxidizer venturi parameters
oxidizer_flow_rate = (thrust / (isp * 9.80665)) * (of_ratio / (1 + of_ratio))
ox_venturi_exit_pressure, ox_venturi_inlet_pressure, ox_required_CdA, ox_venturi_throat_diameter, ox_tube_inner_diameter_oxidizer = calculate_oxidizer_venturi(
    chamber_pressure, dp_injector, line_losses, cavitation_ratio, Cd, oxidizer_flow_rate, oxidizer_fluid, oxidizer_temperature, tube_od_oxidizer, tube_wall_oxidizer
)

# Output results
print("Calculated Venturi Parameters for Methane (Fuel):")
print(f"Venturi Exit Pressure: {venturi_exit_pressure:.2f} psi")
print(f"Venturi Inlet Pressure: {venturi_inlet_pressure:.2f} psi")
print(f"Required CdA: {required_CdA:.6f} m^2")
print(f"Venturi Throat Diameter: {venturi_throat_diameter:.2f} inches")
print(f"Fuel Tube Inner Diameter: {tube_inner_diameter_fuel:.2f} inches")

if venturi_throat_diameter > tube_inner_diameter_fuel:
    print("Warning: The fuel venturi throat diameter exceeds the tube's inner diameter!")

print("\nCalculated Venturi Parameters for LOX (Oxidizer):")
print(f"Venturi Exit Pressure: {ox_venturi_exit_pressure:.2f} psi")
print(f"Venturi Inlet Pressure: {ox_venturi_inlet_pressure:.2f} psi")
print(f"Required CdA: {ox_required_CdA:.6f} m^2")
print(f"Venturi Throat Diameter: {ox_venturi_throat_diameter:.2f} inches")
print(f"Oxidizer Tube Inner Diameter: {ox_tube_inner_diameter_oxidizer:.2f} inches")

if ox_venturi_throat_diameter > ox_tube_inner_diameter_oxidizer:
    print("Warning: The oxidizer venturi throat diameter exceeds the tube's inner diameter!")
