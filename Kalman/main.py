import numpy as np
from Equations import Kalman
import matplotlib.pyplot as plt

A = np.array([[1, 1], [0, 1]])
H = np.array([[1, 0]])
Q = np.array([[0.1, 0], [0, 0.1]])
R = np.array([[0.1]])
P = np.array([[1, 0], [0, 1]])
x0 = np.array([[0], [0]]) #this is gonna have 16 elements instead of 2 lol

kf = Kalman(A, H, Q, R, P, x0)
measurements = [1, 2, 3, 4, 5]

states = []
for z in measurements:
    kf.estimate()
    kf.update(np.array([[z]]))
    states.append(kf.x.flatten()[1])
    print("state:", kf.x.flatten())
plt.plot(measurements, states)
plt.show()
