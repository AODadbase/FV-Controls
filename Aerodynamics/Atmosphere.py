from rocketpy import Environment, SolidMotor, Rocket, Flight
import datetime

#set up presets
#calculating aerodynamic data from preset
def setEnvironment(latitude, longitude, elevation, hour):
    env = Environment(latitude, longitude, elevation)
    tomorrow = datetime.date.today() + datetime.timedelta(1)
    env.set_date((tomorrow.year, tomorrow.month, tomorrow.day, hour)) #hour given in UTC
    env.set_atmospheric_model(type="Forecast", file="GFS")
    URL = "http://weather.uwyo.edu/cgi-bin/sounding?region=samer&TYPE=TEXT%3ALIST&YEAR=2019&MONTH=02&FROM=0500&TO=0512&STNM=83779"
    env.set_atmospheric_model(type="wyoming_sounding", file=URL)
    env.all_info()
