from rocketpy import Environment, SolidMotor, Rocket, Flight
import datetime
class Atmosphere:
    latitude = null
    longitude = null
    elevation = null
    hour = null
    URL = "http://weather.uwyo.edu/cgi-bin/sounding?region=samer&TYPE=TEXT%3ALIST&YEAR=2019&MONTH=02&FROM=0500&TO=0512&STNM=83779"
    

    def __init__ (env, latitude, longitude, elevation, hour):
        env.latitude = latitude
        env.longitude = longitude
        env.elevation = elevation
        env.hour = hour
        today = datetime.date.today()
        env.set_date((today.year, today.month, today.day, hour)) #hour given in UTC
        env.set_atmospheric_model(type="Forecast", file="GFS")
        env.set_atmospheric_model(type="wyoming_sounding", file=URL)


#calculating aerodynamic data from preset
#comment
