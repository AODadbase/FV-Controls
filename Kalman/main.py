import numpy as np
from Equations import Kalman
import matplotlib.pyplot as plt
import random

dT = 0.1
mass = 10
springConst = 4
def initializeValues(variable, ):
    if(variable == "F"):
        return np.array([
            [1,0,0,dT,0,0],
            [0,1,0,0,dT,0],
            [0,0,1,0,0,dT],
            [0,0,0,1,0,0],
            [0,0,0,0,1,0],
            [0,0,0,0,0,1]
        ])
    if(variable == "X"):
        return np.array([[5],[5],[5],[1],[2],[3]]) 
    if(variable == "P"):
        return  np.array([
            [1,0,0,3,0,0],
            [0,1,0,0,3,0],
            [0,0,1,0,0,3],
            [0,0,0,1,0,0],
            [0,0,0,0,1,0],
            [0,0,0,0,0,1]
        ])
    if(variable == "Q"):
        return  0.1 * np.identity(6)
    else:
        return np.zeros(6,6)

def addNoise(value):
    noiseValue = random.randrange(0.9,1.1)
    return value * noiseValue

F = np.array([[1, 1], [0, 1]])
H = np.array([[1, 0]])
Q = np.array([[0.1, 0], [0, 0.1]])
R = np.array([[0.1]])
P = np.array([[1, 0], [0, 1]])
x0 = np.array([[0], [0]]) #this is gonna have 16 elements instead of 2 lol

kf = Kalman(F, H, Q, P, x0, R)
measurements = [1, 2, 3, 4, 5]

#Start with the original one using the other function
#Update both the true value and the kalman one
#Use the measuremnt to do the measuring!

states = []
for z in measurements:
    kf.estimate()
    kf.update(np.array([[z]]))
    states.append(kf.x.flatten()[1])
    print("state:", kf.x.flatten())
plt.plot(measurements, states)
plt.show()


#sliding down a sin wave
#This will be 2D
#For this example, track a oriented frame, 
# State Variables: Position, Velocity, Acceleration, Euler Angle, Change in Euler angles
# vx,vy,ax,ay, theta, theta dot
#Position must be used sparingly, because I'm not sure how much we can use it
#Maybe its better to do a harmonic occilator,


