import numpy as np
import random 
from scipy.optimize import fsolve
import os
#dT = 0.1

defaultMagneticFeild = np.array([1,1,1])



def rotation_matrix_x(theta):
    """Rotation matrix around the X-axis"""
    return np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta), np.cos(theta)]
    ])

def rotation_matrix_y(theta):
    """Rotation matrix around the Y-axis"""
    return np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

def rotation_matrix_z(theta):
    """Rotation matrix around the Z-axis"""
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])




#Convention is RxRyRz * newMagneticVector = ground magneticvector
#So, for the ground, we would multiply by these angles to get what is in the ground
def equations(var, vector):
    xAngle, yAngle, zAngle = var
    rotatedVector = rotation_matrix_x(xAngle) @ rotation_matrix_y(yAngle) @ rotation_matrix_z(zAngle) @ vector - defaultMagneticFeild
    return [rotatedVector[0], rotatedVector[1], rotatedVector[2]]

def solve_for_given_magneticFeild(magneticField):
    func_with_vec = lambda angles : equations(angles, magneticField)
    guess = np.array([np.pi/2, np.pi/2,np.pi/2])
    solution = fsolve(func_with_vec, guess)

    return np.array([solution[0], solution[1], solution[2]])


class Kalman:
    def __init__(self, F, H, Q, P, x0, IT, dIT, mass):
        #F stands for physics, because physics starts with F
        self.F = F 
        #Transforms from nomal basis to input basis
        self.H = H 
        # Q is the uncertainty gain
        self.Q = Q 
        # P is the state uncertainty
        self.P = P 
        # x is the state vector
        self.x = x0 
        # How the inputs affect the stuff
        self.IT = IT
        self.mass = mass
        self.defaultIT = dIT
        self.dT = 0.1

    #Estimates the next step using the change of basis matrix 
    #Inputs are the inputs from the controls
    def estimate(self, inputs):
        self.x = np.dot(self.F, self.x) + np.dot(self.IT, inputs)
        self.P = np.dot(np.dot(self.F, self.P), self.F.T) + self.Q

    #R is the input variance

    def update(self, z_k, R):
        self.updateFMatrix()
        inverted_part = np.dot(np.dot(self.H, self.P), self.H.T) + R
        K = np.dot(np.dot(self.P, self.H.T), np.linalg.inv(inverted_part))
        self.x = self.x + np.dot(K, z_k - np.dot(self.H, self.x))
        self.P = np.dot((np.eye(self.P.shape[0]) - np.dot(K, self.H)), self.P)

    def updateFMatrix(self):
        dT = self.dT
        velo = np.array([self.x[3], self.x[4], self.x[5]])
        DC = 10                                                         #Add in the drag coefficient and AREA PLZ!!!!@@!@@@$@#$#%^&*^$%@#%&*
        area = 10
        DragCoeff = (-0.5 * 1.225 * DC * area * np.linalg.norm(velo))/self.mass
        self.F =np.array([[1,0,0,dT,0,0],
             [0,1,0,0,dT,0],
             [0,0,1,0,0,dT],
             [0,0,0,1 + DragCoeff,0,0],
             [0,0,0,0,1 + DragCoeff,0],
             [0,0,0,0,0,1 + DragCoeff - 9.8/velo[3]]])

    #Need to update the H matrix after every step and the IT matrix
        
    #For this case, the IT matrix will just be adding a force.
    #It's basis will need to be changed
    

    #rotation goes from ground to rocket
    #get rotation matrix angles from magnetometer angles
    def getRotationMatrix(self, magneticField):
        angles = solve_for_given_magneticFeild(magneticField)
        return rotation_matrix_x(angles[0]) @ rotation_matrix_y(angles[1]) @ rotation_matrix_z(angles[2])

    def updateTransitionMatricies(self,measuredField):
      rotationMatrix = self.getRotationMatrix(measuredField)
      
      self.H = np.zeros((6,6))
      for i in range(3):
        for j in range(3):
            if(i == j):
                self.H[i][j] = 1
      for i in range(3):
        for j in range(3):
            self.H[i+3][j+3] = rotationMatrix[i][j]
      self.IT = self.H @ self.defaultIT

    def timeStep(self, inputs, z_k, R, field):
        self.updateTransitionMatricies(field)
        self.estimate(inputs)
        self.update(z_k, R)

def addNoise(value,noise):
    noiseValue = 1 + (random.random() - 0.5)*noise
    return value * noiseValue



class TrueValues:
    def __init__(self,x0, theta0, constants, dt):
        self.x0 = x0
        self.mass = constants["mass"]
        self.Frame = [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
        self.previousPositions = [x0]

        #theta0 stores [theta1, theta2, theat3, w1,w2,w3]
        self.theta = theta0
        self.constants = constants   
        self.dt = dt    
        self.time = 0
        self.thrustT = [] 
        self.thrustS = []
        self.updateFMatrix() 
        self.initializeThrust()
        self.initializePositions()
        
    def initializeThrust(self):
        path = os.path.abspath('.')
        f = open(path + "\\Kalman\\AeroTech_HP-I280DM.eng", "r")
        i = 0
        for line in f:
            if(i > 2):
                vals = line.split()
                self.thrustT.append(float(vals[0]))
                self.thrustS.append(float(vals[1]))
            i += 1
    
    def getThrust(self, time):
        
        if(time > self.thrustT[len(self.thrustT) - 1]):
            return np.array([0,0,0])
        thrustMagnitude = np.interp(time, self.thrustT, self.thrustS)
        return self.Frame[2] * thrustMagnitude

    def initializePositions(self):
        for i in range(5):
            self.x0 = np.dot(self.F, self.x0)
            self.previousPositions.append(self.x0)

    def PropagateForward(self, inputs):
        dT = self.dt
        self.time += dT

        self.updateFMatrix()
        self.x0 = np.dot(self.F, self.x0)
        self.x0[5] -= 9.8 * dT
        thrustA = self.getThrust(self.time)/(self.mass)
        self.x0[3] += thrustA[0] * dT
        self.x0[4] += thrustA[1] * dT
        self.x0[5] += thrustA[2] * dT
        
        velocity1 = (self.x0 - self.previousPositions[-1])/dT
        ThreeVelocity = np.array([velocity1[3], velocity1[4], velocity1[5]])
        Tangent = ThreeVelocity/ np.linalg.norm(ThreeVelocity)
        

        velocity2 = (self.previousPositions[-1] - self.previousPositions[-2])/dT
        ThreeVelocity2 = np.array([velocity2[3], velocity2[4], velocity2[5]])
        acceleration = (ThreeVelocity - ThreeVelocity2)/dT
        BiNormalDirection = np.linalg.cross(ThreeVelocity, acceleration)
        BiNormal = BiNormalDirection / np.linalg.norm(BiNormalDirection)
        Normal = np.linalg.cross(BiNormal, Tangent)

        self.Frame = [Tangent, Normal, BiNormal]
        self.previousPositions.append(self.x0)



        velo = np.array([self.x0[3], self.x0[4], self.x0[5]])
        wPrime = np.array([self.theta[3], self.theta[4]])


        correctiveCoeff = self.constants["Mk"] * np.dot(velo, velo) / np.linalg.norm(wPrime)

        inertias = self.constants["inertias"]
        FinMoments = self.constants["FinMoments"]
        curThe = self.theta
        #Coding in euler's dynamical equations

        w1Dot = (1/inertias[0]) *((inertias[2] - inertias[1]) * curThe[2]*curThe[1] + (correctiveCoeff * curThe[0]) + FinMoments[0])
        w2Dot = (1/inertias[1]) * ((inertias[2] - inertias[0])*curThe[0]*curThe[2] + (correctiveCoeff*curThe[1]) + FinMoments[1])
        w3Dot = (1/inertias[2]) * ((inertias[0] - inertias[1]) *curThe[0]*curThe[1] + FinMoments[2])

        self.theta[3] += w1Dot * dT
        self.theta[4] += w2Dot * dT
        self.theta[5] += w3Dot * dT
        self.theta[0] += self.theta[3] * dT
        self.theta[1] += self.theta[4] * dT
        self.theta[2] += self.theta[5] * dT

        # self.applyWind(dT)
        # self.applyInputs(inputs, dT)

    # def applyWind(self, dt):

    # def applyInputs(self, inputs, dt):

    #FOR WIND, TAKE ELI"S CODE



    #Calculate the right euler angle.
    #Convention is to rotation T -> Z
    # def calculateChangeOfBasisAngles(self, magneticFieldVector):
    #     return solve_for_given_magneticFeild(magneticFieldVector)
    
    def updateFMatrix(self):
        velo = np.array([self.x0[3], self.x0[4], self.x0[5]])
        DC = self.constants["dragCoeff"]                                                         #Add in the drag coefficient and AREA PLZ!!!!@@!@@@$@#$#%^&*^$%@#%&*
        area = self.constants["area"]      
        dT = self.dt
        DragCoeff = (-0.5 * 1.225 * DC * area * np.linalg.norm(velo))/self.mass
        self.F =np.array([[1,0,0,dT,0,0],
                        [0,1,0,0,dT,0],
                        [0,0,1,0,0,dT],
                        [0,0,0,1 + DragCoeff,0,0],
                        [0,0,0,0,1 + DragCoeff,0],
                        [0,0,0,0,0,1 + DragCoeff]])


    def getNoisyState(self):
        positions = self.getPositions()
        velocities = self.getVelocities()
        returnArray = np.array([positions[0],positions[1], positions[2], velocities[0], velocities[1], velocities[2]])
        return returnArray


    def STDToFrameMatrix(self):
        matrix =  np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                matrix[j][i] = self.Frame[i][j]
        return np.linalg.inv(matrix)

    def getFramesMagneticVector(self):
        noisyField = np.array([addNoise(defaultMagneticFeild[0], 0.05),addNoise(defaultMagneticFeild[1], 0.05),addNoise(defaultMagneticFeild[2], 0.05)])
        return self.STDToFrameMatrix() @ noisyField
    
    def getPositions(self):
        return np.array([addNoise(self.x0[0],0.05), addNoise(self.x0[1],0.05), addNoise(self.x0[2], 0.05)])
    
    def getVelocities(self):
        groundVelocities = np.array([addNoise(self.x0[3],0.05), addNoise(self.x0[4],0.05), addNoise(self.x0[5],0.05)])
        return self.STDToFrameMatrix() @ groundVelocities

    def getThetas(self):
        noise = 0.05
        return np.array([addNoise(self.theta[0],noise), addNoise(self.theta[1],noise), addNoise(self.theta[2], noise)])
    
    def getOmegas(self):
        noise = 0.05
        return np.array([addNoise(self.theta[3],noise), addNoise(self.theta[4],noise), addNoise(self.theta[5],noise)])
    
    def getPositionsReal(self):
        return np.array([self.x0[0], self.x0[1], self.x0[2]])
    
    def getVelocitiesReal(self):
        groundVelocities = np.array([self.x0[3], self.x0[4], self.x0[5],0.05])
        return groundVelocities

    def getThetasReal(self):
        return np.array([self.theta[0] % 2 * np.pi, self.theta[1]% 2 * np.pi, self.theta[2]% 2 * np.pi])
    
    def getOmegasReal(self):
        return np.array([self.theta[3], self.theta[4], self.theta[5]])
