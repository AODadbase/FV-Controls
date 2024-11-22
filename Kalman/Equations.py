import numpy as np
import random 
from scipy.optimize import fsolve

dT = 0.1

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

    return np.array(solution[0], solution[1], solution[2])


class Kalman:
    def __init__(self, F, H, Q, P, x0, IT, dIT):
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

        self.defaultIT = dIT
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
    

    #rotation goes from ground to rocket
    #get rotation matrix angles from magnetometer angles
    def getRotationMatrix(self, magneticField):
        angles = solve_for_given_magneticFeild(magneticField)
        return rotation_matrix_x(angles[0]) @ rotation_matrix_y(angles[1]) @ rotation_matrix_z(angles[2])


    def updateTransitionMatricies(self,measuredField):
      rotationMatrix = self.getRotationMatrix(measuredField)
      
      self.H = np.zeros(6,3)
      for i in range(3):
        for j in range(3):
            if(i == j):
                self.H[i][j] = 1
      for i in range(3):
        for j in range(3):
            self.H[i+3][j] = rotationMatrix[i][j]
      self.IT = self.H @ self.defaultIT



    def timeStep(self, inputs, z_k, R, theta1,theta2,theta3):
        self.updateTransitionMatricies(theta1,theta2,theta3)
        self.estimate(inputs)
        self.update(z_k, R)


#TO DO!!!!!!
# We havae ae rocket frenet frame, calculate using velocity and acceleration
#calculate Q, a change of basis matrix to the ground frame since we KNOW what the frame will have a basis in the ground cords
#We know the form of Q because it has to be rotation stuff
#Get the angles and feed the angles into the kalman filter
#The angles will be useful to know where we have to go
#Thats basically it.
#
#
#
#
#
#
#
#

class TrueValues:
    def __init__(self, F,P,x0,IT, Frame):
        self.F = F
        self.P = P
        self.x0 = x0
        self.IT = IT
        self.Frame = [np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])]
        self.previousPositions = [np.array([0,0,0]),np.array([0,0,0]),np.array([0,0,0])]

    def PropagateForward(self, inputs):
        self.x0 = np.dot(self.F, self.x0) + np.dot(self.IT, inputs)

        velocity1 = (self.x0 - self.previousPositions[-1])/dT
        Tangent = velocity1/ np.linalg.norm(velocity1)

        velocity2 = (self.previousPositions[-1] - self.previousPositions[-2])/dT
        acceleration = (velocity1 - velocity2)/dT
        BiNormalDirection = np.linalg.cross(velocity1, acceleration)
        BiNormal = BiNormalDirection / np.linalg.norm(BiNormalDirection)
        Normal = np.linalg.corss(BiNormal, Tangent)

        self.Frame = [Tangent, Normal, BiNormal]
        self.previousPositions.append(self.x0)

    #Calculate the right euler angle.
    #Convention is to rotation T -> Z
    # def calculateChangeOfBasisAngles(self, magneticFieldVector):
    #     return solve_for_given_magneticFeild(magneticFieldVector)
    

    def addNoise(value,noise):
        noiseValue = random.randrange(1-noise,1.1+noise)
        return value * noiseValue

    def getNoisyState(self):
        positions = self.getPositions()
        velocities = self.getVelocities()
        returnArray = np.array([positions[0],positions[1], positions[2], velocities[3], velocities[4], velocities[5]])
        return returnArray


    def STDToFrameMatrix(self):
        matrix =  np.zeros(3,3)
        for i in range(3):
            for j in range(3):
                matrix[j][i] = self.Frame[i][j]
        return np.linalg.inv(matrix)

    def getFramesMagneticVector(self):
        noisyField = np.array(self.addNoise(defaultMagneticFeild[0], 0.05),self.addNoise(defaultMagneticFeild[1], 0.05),self.addNoise(defaultMagneticFeild[2], 0.05))
        return self.STDToFrameMatrix() @ noisyField
    
    def getPositions(self):
        return np.array([self.addNoise(self.x0[0],0.05), self.addNoise(self.x0[1],0.05), self.addNoise(self.x0[2], 0.05)])
    
    def getVelocities(self):
        groundVelocities = np.array([self.addNoise(self.x0[3],0.05), self.addNoise(self.x0[4],0.05), self.addNoise(self.x0[5],0.05)])
        return self.STDToFrameMatrix() @ groundVelocities
