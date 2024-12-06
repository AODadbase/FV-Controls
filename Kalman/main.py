import numpy as np
from Equations import Kalman
from Equations import TrueValues
from eulerEquations import state
import matplotlib.pyplot as plt
import random


angular = [0, 0, 0]
torques = [0, 0, 0]
inertias = [0, 0, 0]
filePath = "C:\\Users\\alber\\Documents\\GitHub\\FV-Controls\\Kalman\\flightVehicleDataV3.csv"

rocket = state(angular, torques, inertias, filePath)
rocket.parseData()
F = np.array([[1,0,0,.1,0,0],
             [0,1,0,0,.1,0],
             [0,0,1,0,0,.1],
             [0,0,0,1,0.01,0],
             [0,0,0,0,1,0.01],
             [0,0,0,0,0,1]])

Fbad = np.array([[1,0,0,.1,0,0],
             [0,1,0,0,.1,0],
             [0,0,1,0,0,.1],
             [0,0,0,1,1,0],
             [0,0,0,0,2,1],
             [0,0,0,0,0,1]])

H = np.zeros((6,3))
Q = np.eye(6)*1.5
IT = np.eye(6)
P = np.eye(6)
x0 = np.array([[0], [0],[0],[0],[1],[500]]) 


kf = Kalman(Fbad, H, Q, P, x0, IT, IT)
truth = TrueValues(F, P, x0, IT)
truth.initializePositions()

inputs = np.array([[0],[0],[0],[0],[0],[0]])
#Start with the original one using the other function
#Update both the true value and the kalman one
#Use the measuremnt to do the measuring!
MUnc = np.eye(6)*100
times = np.arange(100)
statesZ = []
statesX = []
statesY = []
truthVals = []
truthValsX = []
truthValsY = []

for time in range(100):
    truth.PropagateForward(inputs)
    truthState = truth.getNoisyState()
    truthField = truth.getFramesMagneticVector()
    kf.timeStep(inputs, truthState, MUnc, truthField)
    truthVals.append(truthState.flatten()[2])
    truthValsX.append(truthState.flatten()[0])
    truthValsY.append(truthState.flatten()[1])

    statesZ.append(kf.x.flatten()[2])
    statesX.append(kf.x.flatten()[0])
    statesY.append(kf.x.flatten()[1])
    #print("state:", kf.x.flatten())

# plt.plot(times, statesZ, label = 'State')
# plt.plot(times, truthVals, label = 'truth')

# plt.plot(times, statesX, label = 'State')
# plt.plot(times, truthValsX, label = 'truth')


plt.plot(times, statesY, label = 'State')
plt.plot(times, truthValsY, label = 'truth')
plt.legend()
plt.show()


