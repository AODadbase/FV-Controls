import datetime
from rocketpy import Environment, SolidMotor, Rocket, Flight

print("HELLO!")
env = Environment(latitude=41.92298772007185, longitude=-88.06013490408121, elevation=243.43)
tomorrow = datetime.date.today() + datetime.timedelta(days=1)

env.set_date(
    (tomorrow.year, tomorrow.month, tomorrow.day, 12)
)  # Hour given in UTC time
env.set_atmospheric_model(type="Forecast", file="GFS")

#env.info()

ourMOtor = SolidMotor(
    thrust_source="C:\\Users\\alber\\Documents\\GitHub\\FV-Controls\\Kalman\\AeroTech_HP-I280DM.eng",  # Or use a CSV thrust file
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

#ourMOtor.info()


coolRocket = Rocket(
    radius=7.87/200,
    mass=2.259,
    inertia=(0.28, 0.002940, 0.002940),
    power_off_drag=0.560,
    power_on_drag=0.580,
    center_of_mass_without_motor=0.669,
    coordinate_system_orientation="tail_to_nose",
)

coolRocket.add_motor(ourMOtor, position=0.01*(117-86.6))

nose_cone = coolRocket.add_nose(
    length=0.19, kind="lvhaack", position=0.01*(117-0.19)
)

fin_set = coolRocket.add_trapezoidal_fins(
    n=4,
    root_chord=0.203,
    tip_chord=0.0762,
    span=0.0737,
    position=0.01*(117-92.7),
    cant_angle=0,
    sweep_angle=62.8
)

tail = coolRocket.add_tail(
    top_radius=0.0787/2, bottom_radius=0.0572/2, length=0.0381, position=.0381
)



#coolRocket.draw()
#coolRocket.plots.static_margin()

#coolRocket.all_info()

test_flight = Flight(
    rocket=coolRocket, environment=env, rail_length=5.2, inclination=85, heading=0
    )
test_flight.info()