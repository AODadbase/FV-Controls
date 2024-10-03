#Atmosphere
import requests
import time
import matplotlib.pyplot as plt

api_key = '1594000ffa3bac143ac4fbb405b368dc'
location = "Champaign, US"


def get_weather_data(api_key, location):
    url = f"http://api.openweathermap.org/data/2.5/weather?q={location}&appid={api_key}&units=metric"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        weather_data = {
            'location': data['name'],
            'temperature': data['main']['temp'],
            'pressure': data['main']['pressure'],
            'humidity': data['main']['humidity'],
            'wind_speed': data['wind']['speed'],
            'wind_direction': data['wind']['deg'],
            'weather_description': data['weather'][0]['description']
        }
        return weather_data
    else:
        return f"Error: {response.status_code}"
    
weather_data = get_weather_data(api_key,location)

# Constants for the ISA model
T0 = weather_data['temperature']+273.15  # Sea level standard temperature (Kelvin)
P0 = weather_data['pressure']*100  # Sea level standard pressure (Pa)
L = 0.0065   # Temperature lapse rate (K/m)
R = 287.05   # Specific gas constant for dry air (J/(kg·K))
g = 9.80665  # Standard gravity (m/s²)

def isa_atmosphere(altitude):
    """
    Calculate temperature, pressure, and density at a given altitude using the ISA model.
    """
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
    
    return T, P, density

# Example usage
altitudes = list(range(0,5000,500))  # altitudes in meters
temperatures,pressures,densities = [],[],[]
for alt in altitudes:
    T, P, density = isa_atmosphere(alt)
    #print(f"Altitude: {alt} meters -> Temperature: {T:.2f} K, Pressure: {P:.2f} Pa, Density: {density:.4f} kg/m³")
    temperatures.append(T)
    pressures.append(P)
    densities.append(density)


# Plotting the results
fig, ax1 = plt.subplots()

# Temperature plot
ax1.set_xlabel('Altitude (m)')
ax1.set_ylabel('Temperature', color='tab:blue')
ax1.plot(altitudes, temperatures, 'b-', label='Temperature')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# Pressure plot
ax2 = ax1.twinx()
ax2.spines['left'].set_position(('outward', 50))
ax2.set_ylabel('Pressure (Pa)', color='tab:red')
ax2.plot(altitudes, pressures, 'r-', label='Pressure')
ax2.tick_params(axis='y', labelcolor='tab:red')


# Density plot
ax3 = ax2.twinx()  
ax3.spines['right'].set_position(('outward', 30))  # position density axis away from pressure
ax3.set_ylabel('Density (kg/m³)', color='tab:green')
ax3.plot(altitudes, densities, 'g-', label='Density')
ax3.tick_params(axis='y', labelcolor='tab:green')

# Add legends
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
ax3.legend(loc='center right')

# Title and grid
plt.title('Atmospheric Conditions vs Altitude')
plt.grid(True)
plt.tight_layout()
plt.show()

