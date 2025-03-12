import numpy as np
from Equations import Kalman
from Equations import TrueValues
from eulerEquations import state
import matplotlib.pyplot as plt
import random


# dt = 0.1
# F = np.array([[1,0,0,dt,0,0],
#              [0,1,0,0,dt,0],
#              [0,0,1,0,0,dt],
#              [0,0,0,1,0.01,0],
#              [0,0,0,0,1,0.01],
#              [0,0,0,0,0,1]])

# Fbad = np.array([[1,0,0,.1,0,0],
#              [0,1,0,0,.1,0],
#              [0,0,1,0,0,.1],
#              [0,0,0,1,1,0],
#              [0,0,0,0,2,1],
#              [0,0,0,0,0,1]])

# H = np.zeros((6,3))
# Q = np.eye(6)*1.5
# IT = np.eye(6)
# P = np.eye(6)
# x0 = np.array([[0], [0],[0],[0],[1],[500]]) 
# mass =5

# kf = Kalman(Fbad, H, Q, P, x0, IT, IT, mass)
# truth = TrueValues(F, P, x0, IT, mass)
# truth.initializePositions()

# inputs = np.array([[0],[0],[0],[0],[0],[0]])
# #Start with the original one using the other function
# #Update both the true value and the kalman one
# #Use the measuremnt to do the measuring!
# MUnc = np.eye(6)*10
# times = np.arange(100)
# statesZ = []
# statesX = []
# statesY = []
# truthVals = []
# truthValsX = []
# truthValsY = []

# for time in range(100):
#     truth.PropagateForward(inputs)
#     truthState = truth.getNoisyState()
#     truthField = truth.getFramesMagneticVector()
#     kf.timeStep(inputs, truthState, MUnc, truthField)
#     truthVals.append(truthState.flatten()[2])
#     truthValsX.append(truthState.flatten()[0])
#     truthValsY.append(truthState.flatten()[1])

#     statesZ.append(kf.x.flatten()[2])
#     statesX.append(kf.x.flatten()[0])
#     statesY.append(kf.x.flatten()[1])
    #print("state:", kf.x.flatten())

# plt.plot(times, statesZ, label = 'State')
# plt.plot(times, truthVals, label = 'truth')

# plt.plot(times, statesX, label = 'State')
# plt.plot(times, truthValsX, label = 'truth')
print("HELLO")

def initializeParams():
    params = dict()
    params["Mk"] = 5
    params["inertias"] = [2,3,4]
    params["FinMoments"] = [0.01, 0.01, 0.01]
    params["dragCoeff"] = 0.75
    params["area"] = 5
    params["mass"] = 2.5
    return params

constants = initializeParams()

sV = 0.001
x0 = np.array([sV,sV,sV,sV,sV,sV])
thetas = np.array([sV,sV,sV,sV,sV,sV])
simulation = TrueValues(x0, thetas, constants, 0.01)
inputs = np.array([[0],[0],[0],[0],[0],[0]])

thetas = []
pos = []
times = []
for i in range(1000):
    simulation.PropagateForward(inputs)
    thetas.append(simulation.getThetasReal())
    pos.append(simulation.getPositionsReal())
    times.append(simulation.time)

plt.plot(times, pos[:,2], label = 'State')
plt.legend()
plt.show()


