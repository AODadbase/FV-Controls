import sys
import os

# add Folder1 to path
sys.path.append(os.path.abspath("../Control"))


import datetime
import numpy as np
from scipy import linalg

from rocketpy import Environment, SolidMotor, Rocket, Flight
from Dynamics.OurFin import ourFins
from rocketpy.control.controller import _Controller
from Control.ControlSimulation import PhysicsCalc
import matplotlib.pyplot as plt

from Control.ControlSimSympy import Controls


#No constants in yet (I forgot what goes for)
def generateConstants():
    constants = 1
    return constants

#This is the controller that the rocket uses
#Every timestep rollControlFunction is called which will update the aileron angles using Physics Calc
def rollControlFunction(
    time, sampling_rate, state, state_history, observed_variables, finTabs
):
    ABCalculator = PhysicsCalc() #Comes from control Simulation.py
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

def rocketpy_state_to_xhat(state):
    # unpack RocketPy state
    vx, vy, vz   = state[3],  state[4],  state[5]
    e0, e1, e2, e3 = state[6],  state[7],  state[8],  state[9]   # quaternion (scalar-first)
    wx, wy, wz   = state[10], state[11], state[12]

    # your convention is [w1 w2 w3 v1 v2 v3 qw qx qy qz]
    return np.array([wx, wy, wz, vx, vy, vz, e0, e1, e2, e3], dtype=float)

def make_measurement_from_rocketpy(state, sensors=None):
    """
    Returns y_k to feed the observer.
    Choose y to match your C(x,u) output definition.
    Example: measure body angular rates and vertical velocity (4 outputs).
    Replace with ctrl.deriveSensorModel(...) if you already have that.
    """
    wx, wy, wz = state[10], state[11], state[12]
    vz = state[5]

    y = 

    # If you’ve attached RocketPy sensors, you can instead read them:
    # if sensors: y = build_from([s.measurement for s in sensors])
    # Or, if you already implemented: y = ctrl.deriveSensorModel(time, state, sensors)
    return y

def controls_function(
        t: float,
        sampling_rate: float,
        state: list,
        state_history: list,
        observed_variables: list,
        interactive_objects: list,
    ):
    """
    An user-defined function responsible for controlling the simulation.
    This function is expected to take the following arguments, in order:

        1. `time` (float): The current simulation time in seconds.
        2. `sampling_rate` (float): The rate at which the controller
            function is called, measured in Hertz (Hz).
        3. `state` (list): The state vector of the simulation, structured as
            `[x, y, z, vx, vy, vz, e0, e1, e2, e3, wx, wy, wz]`.
        4. `state_history` (list): A record of the rocket's state at each
            step throughout the simulation. The state_history is organized as
            a list of lists, with each sublist containing a state vector. The
            last item in the list always corresponds to the previous state
            vector, providing a chronological sequence of the rocket's
            evolving states.
        5. `observed_variables` (list): A list containing the variables that
            the controller function returns. The return of each controller
            function call is appended to the observed_variables list. The
            initial value in the first step of the simulation of this list is
            provided by the `initial_observed_variables` argument.
        6. `interactive_objects` (list): A list containing the objects that
            the controller function can interact with. The objects are
            listed in the same order as they are provided in the
            `interactive_objects`.
        7. `sensors` (list): A list of sensors that are attached to the
            rocket. The most recent measurements of the sensors are provided
            with the ``sensor.measurement`` attribute. The sensors are
            listed in the same order as they are added to the rocket

    This function will be called during the simulation at the specified
    sampling rate. The function should evaluate and change the interactive
    objects as needed. The function return statement can be used to save
    relevant information in the `observed_variables` list.

    .. note:: The function will be called according to the sampling rate
    specified.
    """
    sampling_rate = 20.0 # Hz
    dt = 1.0 / sampling_rate # seconds
    ## Define gain matrix ##
    Kmax = 100
    Kmin = 17.5

    Ks = np.array([Kmax, Kmin])  # Gain scheduling based on altitude

    ## Define initial conditions ##
    xhat0 = np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0]) # Initial state estimate
    u0 = np.array([0])
    controller = Controls(
        dt=dt,
        Ks=Ks,
        x0=xhat0,
        u0=u0,
    )
    controller.deriveEOM(post_burnout=False)
    controller.deriveEOM(post_burnout=True)
    controller.test_AB(t, rocketpy_state_to_xhat(state), u=u0)
    state_history = controller.states


#Makes a rocket with rocket py
#Values are taken from open rocket
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
    #Remeasure
    ourMOtor = SolidMotor(
    # thrust_source="C:\\Users\\alber\\Documents\\GitHub\\FV-Controls\\Kalman\\AeroTech_HP-I280DM.eng",  # Or use a CSV thrust file
    # thrust_source = '/Users/dsong/Downloads/AeroTech HP-I280DM.eng',
    thrust_source="./Dynamics/AeroTech_HP-I280DM.eng",  # Or use a CSV thrust file
    dry_mass=(0.616 - 0.355),  # kg
    burn_time=1.9,  # Corrected burn time

    dry_inertia=(0.00055, 0.00055, 0.00011),  # kg·m² (approximated)
    nozzle_radius= (10 / 1000), 
    grain_number=5,
    grain_density=18, 
    grain_outer_radius= 7 / 1000,  
    grain_initial_inner_radius=4 / 1000,  
    grain_initial_height= 360 / 5000,  
    grain_separation=0.01,  
    grains_center_of_mass_position=-0.07,  # Estimated
    center_of_dry_mass_position=0.05,  # Estimated
    nozzle_position=-0.3,
    throat_radius= 3.5 / 1000,  
    coordinate_system_orientation="nozzle_to_combustion_chamber",
    )
    coolRocket.add_motor(ourMOtor, position=0.01*(117-86.6))
    nose_cone = coolRocket.add_nose(
        length=0.19, kind="lvhaack", position=0.01*(117-0.19)
    )
    #Boat Tail
    #Verify that it is von karman
    tail = coolRocket.add_tail(
        top_radius=0.0787/2, bottom_radius=0.0572/2, length=0.0381, position=.0381
    )


    #Created in OurFin.py, inherited from rocketpy fins
    ourNewFins = ourFins(
        n=4,
        root_chord=0.203,
        tip_chord=0.0762,
        span=0.0737,
        rocket_radius = 7.87/200,
        cant_angle=0.01,
        sweep_angle=62.8
    )
    #Sampling
    ourController = _Controller(
        interactive_objects= ourNewFins,
        controller_function= rollControlFunction, #Pass our function into rocketpy
        sampling_rate= samplingRate, #How often it runs
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

## To run, run python -m Dynamics.MAURICE2RocketPy