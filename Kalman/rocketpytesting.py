import datetime
import numpy as np
from scipy import linalg

from rocketpy import Environment, SolidMotor, Rocket, Flight
from OurFin import ourFins
from rocketpy.control.controller import _Controller
from eulerEquations import PhysicsCalc
import matplotlib.pyplot as plt
def generateConstants():
    constants = 1
    return constants

def rollControlFunction(
    time, sampling_rate, state, state_history, observed_variables, finTabs
):
    ABCalculator = PhysicsCalc()
    constants = ABCalculator.getConstants(time)

    # constants = generateConstants()
    ourState = np.array([state[10], state[11], state[12], state[3], state[4], state[5]])
    oldAngles = finTabs.aileronAngles
    finTabs.aileronAngles = ABCalculator.getU(time, ourState, constants, oldAngles)
    print("The aileron angles are: "  + str(finTabs.aileronAngles))
    # finTabs.aileronAngles = np.array([0,0,0,0])
    return (
        time,
        finTabs.aileronAngles   
    )




def makeOurRocket(samplingRate):
    coolRocket = Rocket(
        radius=7.87/200,
        mass=2.259,
        inertia=(0.28, 0.002940, 0.002940),
        power_off_drag=0.560,
        power_on_drag=0.580,
        center_of_mass_without_motor=0.669,
        coordinate_system_orientation="tail_to_nose",
    )
    ourMOtor = SolidMotor(
    # thrust_source="C:\\Users\\alber\\Documents\\GitHub\\FV-Controls\\Kalman\\AeroTech_HP-I280DM.eng",  # Or use a CSV thrust file
    thrust_source = '/Users/dsong/Downloads/AeroTech HP-I280DM.eng',
    dry_mass=(0.616 - 0.355),  # kg
    burn_time=1.9,  # Corrected burn time

    dry_inertia=(0.00055, 0.00055, 0.00011),  # kg·m² (approximated)
    nozzle_radius= (10 / 1000),  # 14 mm = 0.014 m
    grain_number=5,
    grain_density=18,  # kg/m³
    grain_outer_radius= 7 / 1000,  # 19 mm = 0.019 m
    grain_initial_inner_radius=4 / 1000,  # 7 mm = 0.007 m
    grain_initial_height= 360 / 5000,  # 120 mm = 0.12 m
    grain_separation=0.01,  # Single grain
    grains_center_of_mass_position=-0.07,  # Estimated
    center_of_dry_mass_position=0.05,  # Estimated
    nozzle_position=-0.3,
    throat_radius= 3.5 / 1000,  # 5.5 mm = 0.0055 m
    coordinate_system_orientation="nozzle_to_combustion_chamber",
    )
    coolRocket.add_motor(ourMOtor, position=0.01*(117-86.6))
    nose_cone = coolRocket.add_nose(
        length=0.19, kind="lvhaack", position=0.01*(117-0.19)
    )
    tail = coolRocket.add_tail(
        top_radius=0.0787/2, bottom_radius=0.0572/2, length=0.0381, position=.0381
    )

    ourNewFins = ourFins(
        n=4,
        root_chord=0.203,
        tip_chord=0.0762,
        span=0.0737,
        rocket_radius = 7.87/200,
        cant_angle=0.01,
        sweep_angle=62.8
    )
    ourController = _Controller(
        interactive_objects= ourNewFins,
        controller_function= rollControlFunction,
        sampling_rate= samplingRate,
        name="KRONOS",
    )

    coolRocket.add_surfaces(ourNewFins, 0.01*(117-92.7))
    coolRocket._add_controllers(ourController)
    return coolRocket, ourController





#Initializing the rocket simulation
print("HELLO!")
env = Environment(latitude=41.92298772007185, longitude=-88.06013490408121, elevation=243.43)
tomorrow = datetime.date.today() + datetime.timedelta(days=1)
env.set_date((tomorrow.year, tomorrow.month, tomorrow.day, 12))  
env.set_atmospheric_model(type="Forecast", file="GFS")

coolRocket, ourController = makeOurRocket(100)






#coolRocket.draw()
#coolRocket.plots.static_margin()
#coolRocket.all_info()




#Test Flight
test_flight = Flight(
    rocket=coolRocket, environment=env, rail_length=5.2, inclination=85, heading=0
    )
test_flight.info()
test_flight.plots.angular_kinematics_data()

test_flight.export_data(
    "testing.csv",
    "w1",
    "w2",
    "w3",
    "alpha1",
    "alpha2",
    "alpha3",

)

#https://github.com/RocketPy-Team/RocketPy/blob/master/rocketpy/control/controller.py
# Rocket py controller class, so that you can control stuff

#https://docs.rocketpy.org/en/latest/reference/classes/aero_surfaces/Fins.html
#Fins.roll_parameters (list) – List containing the roll moment lift coefficient, 
#the roll moment damping coefficient and the cant angle in radians.

#can create a custom fin class, and change the moments, I think then I can make a controller object


##I COMMENTED OUT THE CONTROLLER AND COMMENTED OUT THE PHYSICS CHANGES IN OURFIN