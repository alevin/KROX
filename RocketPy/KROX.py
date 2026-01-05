from math import exp 

from rocketpy import units
from rocketpy import Fluid 
from rocketpy import LiquidMotor
from rocketpy import CylindricalTank, MassFlowRateBasedTank

OXmDot_engine = 1 #kg/sec
KRmDot_engine = 0.5 #kg/sec
burnDuration_engine = 10 #seconds with OX running out first
thurst_engine = units.convert_units(1200, "lb", "kg")
dryweight_engine = units.convert_units(2, "lb", "kg")
nozzleOutletRadius_engine = 0.075 #Radius of motor nozzle outlet in meters.
centerOfDryMassPositionRelativeToNozzle_engine = units.convert_units(8, "in", "m") #exit of nozzle is at origin.

#there is no pressurant yet

# LOX tank parameters
tankRadiusOX = units.convert_units(5.562/2, "in", "m")
tankLengthOX = units.convert_units(32, "in", "m")

# kerosene tank parameters
tankRadiusKR = units.convert_units(5.562/2, "in", "m")
tankLengthKR = units.convert_units(32, "in", "m")




#creating tank geometry objects. standard ones here, not detailed. 
OXtank_geometry = CylindricalTank(radius=0.1, height=2.0, spherical_caps=True)
KRtank_geometry = CylindricalTank(radius=0.1, height=2.0, spherical_caps=True)

#defining our fluids, kerosene, lox, air. These values are approximated from engineering toolbox
#this uh, more accurate please LOL
OX_liquid = Fluid(name="Liquid Oxygen", density=1000)  #kg/m^3 somehwere near -100degC, somewhere near 50bara 
OX_vapour = Fluid(name="Vapour Oxygen", density=26)    #kg/m^3 somehwere near -100degC, somewhere near 10bara 

KR_liquid = Fluid(name="Liquid Kerosene", density=800) #kg/m^3 range was given as 775/840 o.o
AIR_vapour = Fluid(name="Air", density=10)              #kg/m^3 at 20C? bleh

#last i checked, the engine is eating at a constant m-dot, so mass flow tanks it will be

OXtank_fillMass = OX_liquid.density * OXtank_geometry.total_volume #assuming 100% flled initial tank
KRtank_fillMass = KR_liquid.density * KRtank_geometry.total_volume #assuming 100% flled initial tank

KRtank_flameballtime = (KRtank_fillMass - burnDuration_engine * KRmDot_engine) / KRmDot_engine #firing is OX limited. this figuring out leftover time the tank is dumping fuel
print("going")


# Define fluids
oxidizer_liq = Fluid(name="N2O_l", density=1220)
oxidizer_gas = Fluid(name="N2O_g", density=1.9277)
fuel_liq = Fluid(name="ethanol_l", density=789)
fuel_gas = Fluid(name="ethanol_g", density=1.59)

# Define tanks geometry
tanks_shape = CylindricalTank(radius = 0.1, height = 1.2, spherical_caps = True)

# Define tanks
oxidizer_tank = MassFlowRateBasedTank(
    name="oxidizer tank",
    geometry=tanks_shape,
    flux_time=5,
    initial_liquid_mass=32,
    initial_gas_mass=0.01,
    liquid_mass_flow_rate_in=0,
    liquid_mass_flow_rate_out=lambda t: 32 / 3 * exp(-0.25 * t),
    gas_mass_flow_rate_in=0,
    gas_mass_flow_rate_out=0,
    liquid=oxidizer_liq,
    gas=oxidizer_gas,
)

fuel_tank = MassFlowRateBasedTank(
    name="fuel tank",
    geometry=tanks_shape,
    flux_time=5,
    initial_liquid_mass=21,
    initial_gas_mass=0.01,
    liquid_mass_flow_rate_in=0,
    liquid_mass_flow_rate_out=lambda t: 21 / 3 * exp(-0.25 * t),
    gas_mass_flow_rate_in=0,
    gas_mass_flow_rate_out=lambda t: 0.01 / 3 * exp(-0.25 * t),
    liquid=fuel_liq,
    gas=fuel_gas,
)











"""
OXtank = MassFlowRateBasedTank(
    name="Liquid Oxygen Tank",
    geometry=OXtank_geometry,    #defined earlier
    flux_time=burnDuration_engine-1,          #engine burn time basically.
    liquid=OX_liquid,            #defined earlier
    gas=OX_vapour,               #hey so, is it just ox vapour? wouldnt it be both? someone run to mixed gas properties
    initial_liquid_mass=OXtank_fillMass,    #earlier
    initial_gas_mass=0,       #smol amount 
    liquid_mass_flow_rate_in=0,  #should be
    liquid_mass_flow_rate_out=OXmDot_engine,
    gas_mass_flow_rate_in=0.15,     #well. how much...?
    gas_mass_flow_rate_out=0.05,    #well. none zero due to boil off. 
    discretize=100,              # i did not understand this
)

KRtank = MassFlowRateBasedTank(
    name="Kerosene Tank",
    geometry=KRtank_geometry,    #defined earlier
    flux_time=burnDuration_engine+KRtank_flameballtime,     #how long does tank need to be worried about
    liquid=KR_liquid,            #defined earlier
    gas=AIR_vapour,              #all air
    initial_liquid_mass=KRtank_fillMass,    #earlier
    initial_gas_mass=0,       #smol amount 
    liquid_mass_flow_rate_in=0,  #should be
    liquid_mass_flow_rate_out=KRmDot_engine,
    gas_mass_flow_rate_in=0.3,     #well. how much...?
    gas_mass_flow_rate_out=0,    #basically 0. 
    discretize=100,              # i did not understand this
)

KROXengine = LiquidMotor(
    coordinate_system_orientation="nozzle_to_combustion_chamber", #positive torwards to combustion chamber
    center_of_dry_mass_position=centerOfDryMassPositionRelativeToNozzle_engine,
    nozzle_position=0, #nozzle exit area to coord sys
    thrust_source= thurst_engine,
    dry_mass=dryweight_engine,
    dry_inertia=(0.125, 0.125, 0.002),           #uhm, need this.
    nozzle_radius=nozzleOutletRadius_engine,
    burn_time=burnDuration_engine,
)
KROXengine.add_tank(tank=OXtank, position=1.0)
KROXengine.add_tank(tank=KRtank, position=2.5)

KROXengine.info()
"""