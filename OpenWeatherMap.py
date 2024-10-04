#Atmosphere
import requests
import time
from pyhigh import get_elevation


class atmosphere:
    def __init__(self,api_key,location):
        self.api_key = api_key
        self.location = location
        self.temperatures = []
        self.pressures = []
        self.densities = []
        self.altitudes = []


    def get_weather_data(self):
        url = f"http://api.openweathermap.org/data/2.5/weather?q={self.location}&appid={self.api_key}&units=metric"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            #print(data)
            self.lon = data['coord']['lon']
            self.lat = data['coord']['lat']
            self.weather_data = {
                'location': data['name'],
                'temperature': data['main']['temp'],
                'pressure': data['main']['pressure'],
                'humidity': data['main']['humidity'],
                'wind_speed': data['wind']['speed'],
                'wind_direction': data['wind']['deg'],
                'weather_description': data['weather'][0]['description']
            }
            return self.weather_data
        else:
            return f"Error: {response.status_code}"
        
    
    def isa_atmosphere(self,max_altitude):
        # Constants for the ISA model
        min_altitude = int(round(get_elevation(self.lat,self.lon),-1))
        #clear_cache() #optional, saves computer space however makes program take longer
        T0 = self.weather_data['temperature']+273.15  # Sea level standard temperature (Kelvin)
        P0 = self.weather_data['pressure']*100  # Sea level standard pressure (Pa)
        L = 0.0065   # Temperature lapse rate (K/m)
        R = 287.05   # Specific gas constant for dry air (J/(kg·K))
        g = 9.80665  # Standard gravity (m/s²)
        for altitude in range(min_altitude,max_altitude+10,10):
            if altitude < 11000:  # Troposphere (up to 11 km)
                T = T0 - L * altitude
                P = P0 * (T / T0) ** (-g / (R * L))
            else:
                T = None
                P = None         
            if T and P:
                # Calculate air density using the ideal gas law: ρ = P / (R * T)
                density = P / (R * T)
            else:
                density = None
            self.altitudes.append(altitude)
            self.temperatures.append(round(T,3))
            self.pressures.append(round(P,3))
            self.densities.append(round(density,3))


    def print_data(self):
        print(self.altitudes[0])
        print(self.pressures[0])
        print(self.densities[0])
        print(self.temperatures[0])


    def graph(self,p=True,t=True,d=True):
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots()
        if p: #pressure
            ax1.set_xlabel('Altitude (m)')
            ax1.set_ylabel('Pressure (Pa)', color='tab:blue')
            ax1.plot(self.altitudes, self.pressures, 'b-', label='Pressure')
            ax1.tick_params(axis='y', labelcolor='tab:blue')
            ax1.legend(loc='upper left')
        if t: #temperature
            ax2 = ax1.twinx()
            ax2.spines['left'].set_position(('outward', 50))
            ax2.set_ylabel('Temperature (C)', color='tab:red')
            ax2.plot(self.altitudes, self.temperatures, 'r-', label='Pressure')
            ax2.tick_params(axis='y', labelcolor='tab:red')
            ax2.legend(loc='upper right')
        if d: #density
            ax3 = ax2.twinx()  
            ax3.spines['right'].set_position(('outward', 30))  # position density axis away from pressure
            ax3.set_ylabel('Density (kg/m³)', color='tab:green')
            ax3.plot(self.altitudes, self.densities, 'g-', label='Density')
            ax3.tick_params(axis='y', labelcolor='tab:green')
            ax3.legend(loc='center right')
        plt.title('Atmospheric Conditions vs Altitude')
        plt.grid(True)
        plt.tight_layout()
        plt.show()


    def returndata(self,p=True,t=True,d=True):
        output = [self.altitudes]
        if p:
            output.append(self.pressures)
        if t:
            output.append(self.temperatures)
        if d:
            output.append(self.densities)
        return output



#main process
def getdata(atmosphere,max_altitude):
    atmosphere.get_weather_data()
    atmosphere.isa_atmosphere(max_altitude)
    #atmosphere.print_data()
    print(atmosphere.returndata())
    atmosphere.graph()
    return(atmosphere.returndata())
    

api_key = '1594000ffa3bac143ac4fbb405b368dc'
location = "Champaign, US"
champaignweather = atmosphere(api_key,location)
getdata(champaignweather,10000)

