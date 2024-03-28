from rocketpy import Environment, SolidMotor, Rocket, Flight
#set up presets
#calculating aerodynamic data from preset
Environment.latitude (float) # Launch site latitude.
Environment.longitude (float) # Launch site longitude.
Environment.elevation (float) # Launch site elevation.
Environment.elev_array (array) # Two-dimensional Array containing the elevation information
Environment.pressure (Function) # Air pressure in Pa as a function of altitude. Can be accessed as regular array, or called as a Function. 

wyoming_sounding #sets pressure, temperature, wind-u and wind-v profiles and surface elevation obtained from an upper air sounding given by the file parameter through an URL. This URL should point to a data webpage given by selecting plot type as text: list, a station and a time at weather.uwyo. 
