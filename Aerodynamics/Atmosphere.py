from rocketpy import Environment, SolidMotor, Rocket, Flight
import datetime
class Atmosphere:
    
    URL = "https://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=X^X^X^X^&MONTH=X^X^&FROM=1912&TO=1912&STNM=X^X^X^X^X^"
    

    def __init__ (self, env, latitude, longitude, elevation, hour, date):
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation
        self.hour = hour
        today = datetime.date.today()
        self.date = Environment.set_date((today.year, today.month, today.day, hour)) #hour given in UTC
        self.env = Environment(self.date, self.latitude, self.longitude, self.elevation)
        self.set_atmospheric_model(type="Forecast", file="GFS")
        self.set_atmospheric_model(type="wyoming_sounding", file=URL)

Bogus = Atmosphere(2, 3, 3, 3)
print(Bogus.URL)


#calculating aerodynamic data from preset
