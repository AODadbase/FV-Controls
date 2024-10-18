import numpy as np
import random 
dT = 0.1
class Kalman:
    def __init__(self, F, H, Q, P, x0, IT):
        #F stands for physics, because physics starts with F
        self.F = F 
        #Transforms between the two basis
        self.H = H 
        # Q is the uncertainty gain
        self.Q = Q 
        # P is the state uncertainty
        self.P = P 
        # x is the state vector
        self.x = x0 
        # How the inputs affect the stuff
        self.IT = IT
        self.OldPoistions = self.x

    #Estimates the next step using the change of basis matrix 
    #Inputs are the inputs from the controls
    def estimate(self, inputs):
        self.x = np.dot(self.F, self.x) + np.dot(self.IT, inputs)
        self.P = np.dot(np.dot(self.F, self.P), self.F.T) + self.Q

    #R is the input variance

    def update(self, z_k, R):
        self.OldPoistions = self.x
        inverted_part = np.dot(np.dot(self.H, self.P), self.H.T) + R
        K = np.dot(np.dot(self.P, self.H.T), np.linalg.inv(inverted_part))
        self.x = self.x + np.dot(K, z_k - np.dot(self.H, self.x))
        self.P = np.dot((np.eye(self.P.shape[0]) - np.dot(K, self.H)), self.P)


    #Need to update the H matrix after every step and the IT matrix
        
    #For this case, the IT matrix will just be adding a force.
    #It's basis will need to be changed

    #get rotation matrix angles from magnetometer angles
    def updateTransitionMatricies(self,theta1,theta2,theta3):
        c1 = np.cos(theta1)
        c2 = np.cos(theta2)
        c3 = np.cos(theta3)
        s1 = np.sin(theta1)
        s2 = np.sin(theta2)
        s3 = np.sin(theta3)
        xRotation = np.array([
            [1,0,0,0,0,0],
            [0,1,0,0,0,0],
            [0,0,1,0,0,0],
            [0,0,0,1,0,0],
            [0,0,0,0,c1,-1 *s1],
            [0,0,0,0,s1,c1],
        ])
        yRotation = np.array([
            [1,0,0,0,0,0],
            [0,1,0,0,0,0],
            [0,0,1,0,0,0],
            [0,0,0,c2,0,s2],
            [0,0,0,0,1,0],
            [0,0,0,-1* s2,0,c1],
        ])
        zRotation = np.array([
            [1,0,0,0,0,0],
            [0,1,0,0,0,0],
            [0,0,1,0,0,0],
            [0,0,0,c3,-1 * s3,0],
            [0,0,0,s3,c3,0],
            [0,0,0,0,0,1],
        ])
        self.H = np.dot(np.dot(zRotation, yRotation),xRotation)
        self.x0[8] = (self.x0[7] - self.OldPoistions[7])/dT
        self.IT = np.array([
            [0,0,0],
            [0,0,0],
            [0,0,0],
            [c2*c3,s1*s2*c3 - c1*s3,c1*s3 + c1*s2*c3],
            [c2*s3,c3*c1 + s1*s2*s3,-c3*c1 + c1*s2*s3],
            [-s2,c1*c2,c1*c2]
        ])

    def timeStep(self, inputs, z_k, R, theta1,theta2,theta3):
        self.updateTransitionMatricies(theta1,theta2,theta3)
        self.estimate(inputs)
        self.update(z_k, R)


class TrueValues:
    def __init__(self, F,P,x0,IT):
        self.F = F
        self.P = P
        self.x0 = x0
        self.IT = IT

    def PropagateForward(self, inputs):
        oldX = self.x0
        self.x = np.dot(self.F, self.x) + np.dot(self.IT, inputs)
        if(self.x0[0] == 0) :
            self.x0[7] = 0
        if(self.x0[0] < 0):
            self.x0[7] = 
        if(self.x0[0] > 0):
            self.x0[7] = 

    #Calculate the right euler angle.

    def addNoise(value,noise):
        noiseValue = random.randrange(1-noise,1.1+noise)
        return value * noiseValue


    def getNoisyValues(self,value,noise):
        if(value == "X"):
            a = []
            for x in self.x0:
                a.append(self.addNoise(x,noise))
            return np.array(a)

        if(value == "P"):
            val = 0.05
            return np.identity(6) * noise