"""KROX Rocket Engine Simulation using RocketPy.

This module simulates a liquid oxygen (LOX) and kerosene engine with
pressure-fed propellant tanks.
"""
from math import exp
from typing import Final
from dataclasses import dataclass

from rocketpy import units
from rocketpy import Fluid
from rocketpy import LiquidMotor
from rocketpy import CylindricalTank, MassFlowRateBasedTank


@dataclass
class EngineConfig:
    """Engine configuration parameters."""
    ox_mdot: float = 1.0  # kg/s - Oxidizer mass flow rate
    kr_mdot: float = 0.5  # kg/s - Kerosene mass flow rate
    burn_duration: float = 10.0  # seconds (OX limited)
    thrust_lbf: float = 1200.0  # lbf
    dry_weight_lb: float = 2.0  # lb
    nozzle_outlet_radius: float = 0.075  # meters
    center_of_dry_mass_in: float = 8.0  # inches from nozzle exit


# Initialize engine configuration
config = EngineConfig()

# Engine parameters - mass flow rates used for propellant consumption and pressurant gas
ox_mdot_engine: Final[float] = config.ox_mdot
kr_mdot_engine: Final[float] = config.kr_mdot

# Engine performance parameters
burn_duration_engine: Final[float] = config.burn_duration
thrust_engine: Final[float] = config.thrust_lbf  # RocketPy uses lbf natively
dry_weight_engine: Final[float] = units.convert_units(config.dry_weight_lb, "lb", "kg")
nozzle_outlet_radius_engine: Final[float] = config.nozzle_outlet_radius
center_of_dry_mass_position_engine: Final[float] = units.convert_units(
    config.center_of_dry_mass_in, "in", "m"
)

# Tank geometry parameters
# LOX tank dimensions
tank_radius_ox: Final[float] = units.convert_units(5.562 / 2, "in", "m")
tank_length_ox: Final[float] = units.convert_units(32, "in", "m")

# Kerosene tank dimensions
tank_radius_kr: Final[float] = units.convert_units(5.562 / 2, "in", "m")
tank_length_kr: Final[float] = units.convert_units(32, "in", "m")

# Create tank geometry objects with spherical end caps
ox_tank_geometry = CylindricalTank(
    radius=tank_radius_ox,
    height=tank_length_ox,
    spherical_caps=True
)
kr_tank_geometry = CylindricalTank(
    radius=tank_radius_kr,
    height=tank_length_kr,
    spherical_caps=True
)

# Define fluids - values approximated from Engineering Toolbox
# TODO: Verify densities at actual operating conditions and update
ox_liquid = Fluid(
    name="Liquid Oxygen",
    density=1141  # kg/m³ at -183°C, 1 atm (actual: verify at operating pressure ~50 bar)
)
ox_vapor = Fluid(
    name="Gaseous Oxygen",
    density=26  # kg/m³ (estimate at -100°C, ~10 bar)
)

kr_liquid = Fluid(
    name="Liquid Kerosene",
    density=810  # kg/m³ at 15°C (typical range: 775-840 kg/m³)
)
air_vapor = Fluid(
    name="Air",
    density=10  # kg/m³ (pressurized ullage gas, ~8 bar)
)

# Calculate initial propellant masses
# Use 98% volume for liquid, 2% volume for gas ullage
ox_tank_liquid_volume: float = ox_tank_geometry.total_volume * 0.98
ox_tank_gas_volume: float = ox_tank_geometry.total_volume * 0.02
ox_tank_fill_mass: float = ox_liquid.density * ox_tank_liquid_volume
ox_tank_initial_gas_mass: float = air_vapor.density * ox_tank_gas_volume

kr_tank_liquid_volume: float = kr_tank_geometry.total_volume * 0.98
kr_tank_gas_volume: float = kr_tank_geometry.total_volume * 0.02
kr_tank_fill_mass: float = kr_liquid.density * kr_tank_liquid_volume
kr_tank_initial_gas_mass: float = air_vapor.density * kr_tank_gas_volume

# Validate propellant masses
if ox_tank_fill_mass <= 0 or kr_tank_fill_mass <= 0:
    raise ValueError("Tank fill masses must be positive. Check tank geometry and fluid density.")

# Calculate additional burn time for kerosene after OX depletion ("flameball" time)
# Engine burn is OX-limited, so KR will continue flowing briefly
kr_tank_flameball_time: float = (
    kr_tank_fill_mass - burn_duration_engine * kr_mdot_engine
) / kr_mdot_engine


def calculate_pressurant_flow_rate(
    liquid_mdot: float,
    liquid_density: float,
    gas_density: float
) -> float:
    """Calculate pressurant gas mass flow rate to displace liquid volume.
    
    As liquid flows out, pressurant gas must flow in to maintain pressure
    and fill the displaced volume.
    
    Args:
        liquid_mdot: Liquid mass flow rate (kg/s)
        liquid_density: Liquid density (kg/m³)
        gas_density: Pressurant gas density (kg/m³)
    
    Returns:
        Gas mass flow rate (kg/s)
    """
    volume_flow_rate = liquid_mdot / liquid_density
    return volume_flow_rate * gas_density


# Calculate pressurant gas flow rates for each tank
kr_gas_mdot: float = calculate_pressurant_flow_rate(
    kr_mdot_engine, kr_liquid.density, air_vapor.density
)
ox_gas_mdot: float = calculate_pressurant_flow_rate(
    ox_mdot_engine, ox_liquid.density, air_vapor.density
)



# Create oxidizer tank
ox_tank = MassFlowRateBasedTank(
    name="Liquid Oxygen Tank",
    geometry=ox_tank_geometry,
    flux_time=burn_duration_engine - 1,  # Simulation time window
    liquid=ox_liquid,
    gas=air_vapor,  # TODO: Consider using ox_vapor for ullage gas instead of air
    initial_liquid_mass=ox_tank_fill_mass,
    initial_gas_mass=ox_tank_initial_gas_mass,
    liquid_mass_flow_rate_in=0,  # No refilling during burn
    liquid_mass_flow_rate_out=ox_mdot_engine,
    gas_mass_flow_rate_in=0,  # RocketPy handles ullage gas expansion internally
    gas_mass_flow_rate_out=0,  # Neglecting boiloff losses
    discretize=100,  # Number of time steps for simulation accuracy
)

# Create kerosene tank
kr_tank = MassFlowRateBasedTank(
    name="Kerosene Tank",
    geometry=kr_tank_geometry,
    flux_time=burn_duration_engine + kr_tank_flameball_time,  # Total tank operation time
    liquid=kr_liquid,
    gas=air_vapor,  # Pressurant gas in ullage space
    initial_liquid_mass=kr_tank_fill_mass,
    initial_gas_mass=kr_tank_initial_gas_mass,
    liquid_mass_flow_rate_in=0,  # No refilling during burn
    liquid_mass_flow_rate_out=kr_mdot_engine,
    gas_mass_flow_rate_in=0,  # RocketPy handles ullage gas expansion internally
    gas_mass_flow_rate_out=0,  # Negligible gas losses
    discretize=100,  # Number of time steps for simulation accuracy
)

# Create liquid engine assembly
krox_engine = LiquidMotor(
    coordinate_system_orientation="nozzle_to_combustion_chamber",
    center_of_dry_mass_position=center_of_dry_mass_position_engine,
    nozzle_position=0,  # Nozzle exit at coordinate system origin
    thrust_source=thrust_engine,
    dry_mass=dry_weight_engine,
    dry_inertia=(0.125, 0.125, 0.002),  # TODO: Calculate actual inertia values
    nozzle_radius=nozzle_outlet_radius_engine,
    burn_time=burn_duration_engine,
)

# Add propellant tanks to engine
# Position values are distances from nozzle exit in meters
krox_engine.add_tank(tank=ox_tank, position=1.0)
krox_engine.add_tank(tank=kr_tank, position=2.5)


if __name__ == "__main__":
    # Display engine configuration and performance data
    krox_engine.info()