import numpy as np
import csv
from scipy.optimize import approx_fprime
import sympy as sp
from scipy import linalg

class PhysicsCalc:
    def __init__(self):
        self.something = 1
    #Returns constants, some of this isn't real values
    #Fin moments aren't real values
    #Drag coffecients and inertias are from openrocket, as a function of time

    #CHECK ALL VALUES!!!!!!!!!!!!!!!!!!
    
    def getConstants(self, time):
        constants = dict()
        Inertia = np.array([3,3,0.004])
        constants["Fin Moments"] = np.array([0.003,0.003,0.003]) #Fin misalignment from open rocket
        constants["correctiveConstants"] = 1

        if(time > 1.97):
            constants["Mass"] = 2.991
            Inertia[0] = 0.345
            Inertia[1] = 0.345
            constants["DragCoeff"] = 0.668 -0.0144* time + 0.00121 * time * time - 0.0000353 * time * time * time

        else:
            constants["Mass"] = 2.375
            longI  = 0.32 + (0.32 - 0.28) * (-1 * time) * 0.5
            Inertia[0] = longI
            Inertia[1] = longI
            constants["DragCoeff"] = 0.611 - 0.0155* time + 0.0337 * time * time - 0.00823 * time * time * time



        constants["Inertias"] = Inertia
        return constants

    #This is a helper fucntion for calculating the A matrix
    # Takes in x and returns xdot in xdot = Ax
    def aLambdaFunc(self, x, FinMoments,Mk, Dk , mass ,inertias):
        #x = (w1, w2, w3, v1, v2, v3)
        velo = np.array([x[3], x[4], x[5]])
        wPrime = np.array([x[0], x[1]])


        correctiveCoeff = Mk * np.dot(velo, velo) / np.linalg.norm(wPrime)
        vmag = np.linalg.norm(velo)

        w1Dot = (1/inertias[0]) *((inertias[2] - inertias[1]) * x[2]*x[1] + (correctiveCoeff * x[0]) + FinMoments[0])
        w2Dot = (1/inertias[1]) * ((inertias[2] - inertias[0])*x[0]*x[2] + (correctiveCoeff*x[1]) + FinMoments[1])
        w3Dot = (1/inertias[2]) * ((inertias[0] - inertias[1]) *x[0]*x[1] + FinMoments[2])
        v1Dot = -(Dk * vmag * x[3])/mass
        v2Dot = -(Dk * vmag * x[4])/mass
        v3Dot = -(Dk * vmag * x[5])/mass - 9.8

        return np.array([w1Dot, w2Dot,w3Dot, v1Dot, v2Dot, v3Dot ])
    
    #Uses the aLambdaFunc to calculate A using a lambda function
    def calculateANew(self, constants, state):
        FinMoments = constants["Fin Moments"]
        Mk = constants["correctiveConstants"]
        Dk = constants["DragCoeff"]
        mass = constants["Mass"]
        inertias = constants["Inertias"]
        
        physicsWrapped = lambda x : self.aLambdaFunc(x,FinMoments, Mk, Dk, mass, inertias)
        epsilon = 1e-6
        #approx_fprim calculates the jacobian
        A = approx_fprime(state, physicsWrapped, epsilon)
        return A

    #Calculates the moment. Imput data from CFD
    def calculateMoments(self, velocities, alphas):
        vScale = np.sqrt(np.dot(velocities,velocities))/210
        # Mx = 0.5 *( alphas[0] + alphas[2]) * vScale
        # My =  0.5*( alphas[1] + alphas[3]) * vScale
        
        # Mz = 0.02 * (alphas[1] + alphas[3] - alphas[2] - alphas [0]) * vScale
        Mx = 0.5 *( alphas[0]) * vScale * vScale
        My = 0
        Mz = 1 * (alphas [0]) * vScale * vScale

        return np.array([Mx, My, Mz])

    def calculateBFunctional(self, alphas, x, inertias):
        moments = self.calculateMoments(np.array([x[3], x[4], x[5]]), alphas)
        w1Dot = moments[0]/inertias[0]
        w2Dot = moments[1]/inertias[1]
        w3Dot = moments[2]/inertias[2]
        return (w1Dot, w2Dot, w3Dot,w1Dot * 1e-8,w2Dot * 1e-8,w3Dot * 1e-8)
    
    def calculateBNew(self, alphas, constants, state):
        inputsWrapped = lambda alpha : self.calculateBFunctional(alpha, state, constants["Inertias"])
        epsilon = 1e-6
        #approx_fprim calculates the jacobian
        B = approx_fprime(alphas, inputsWrapped, epsilon)
        return B
        
    def getU(self,time, state, constants, oldAngles):
        # print(state)
        # if(time < 1):
        #     return [0,0,0,0]
            # ourState = np.array([state[3], state[4], state[5], state[10], state[11], state[12]])
        AMatrix = self.calculateANew(constants, state)
        BMatrix = self.calculateBNew(oldAngles, constants,state)

        Qmatrix = np.eye(6)
        Qmatrix[3][3] = 1e-6
        Qmatrix[4][4] = 1e-6
        Qmatrix[5][5] = 1e-6
        Qmatrix[0][0] = 1e-6
        Qmatrix[1][1] = 1e-6
        Qmatrix[2][2] = 100

        Rmatrix = 1/(0.13962634016**2) *  np.eye(4)

        # print(AMatrix)
        # print(BMatrix)
        # try:
        #     P = linalg.solve_continuous_are(AMatrix, BMatrix, Qmatrix, Rmatrix) 
        #     K = linalg.inv(Rmatrix) @ BMatrix.T @ P
        #     print(K)
        # except Exception as e:
        #     print("❌ CARE solve failed:", e)
        #     return [0, 0, 0, 0]

        K = np.array([[0,0,5,0,0,0],
                     [0,0,0,0,0,0],
                     [0,0,0,0,0,0],
                     [0,0,0,0,0,0]])

        #u = -kx
        # for i in range(3):
        #     for j in range(4):
        #         K[j][i + 3] = 0
        u = -1 * K @ state

        # print(K)
        # print("Time is " + str(time))
        # print(u)
        i = 0
        for val in u:
            if(np.abs(val) > 19 * np.pi/180):
                u[i] = np.sign(val) * 10 * np.pi/180
            i = i + 1
        return u





def symPyJacobian():
    w1, w2, w3, v1, v2, v3 = sp.symbols('w1 w2 w3 v1 v2 v3')
    fMx, fMy, fMz = sp.symbols("fMx fMy fMz")
    Mk, Dk, mass = sp.symbols("Mk Dk mass") 
    I1, I2, I3 = sp.symbols("I1 I2 I3")


    correctiveCoeff = Mk * (v1 * v1 + v2 * v2 + v3 * v3) / (sp.sqrt(w1 * w1 + w2 * w2))
    vmag = sp.sqrt(v1 * v1 + v2 * v2 + v3 * v3)

    w1Dot = (1/I1) *((I3 - I2) * w3*w2 + (correctiveCoeff * w1) + fMx)
    w2Dot = (1/I2) * ((I3 - I1)*w1*w3 + (correctiveCoeff*w2) + fMy)
    w3Dot = (1/I3) * ((I1 - I2) *w1*w2 + fMz)
    v1Dot = (Dk * vmag * v1)/mass
    v2Dot = (Dk * vmag * v2)/mass
    v3Dot = (Dk * vmag * v3)/mass - 9.8

    F = sp.Matrix([w1Dot, w2Dot, w3Dot, v1Dot, v2Dot, v3Dot])
    J = F.jacobian([w1, w2, w3, v1, v2, v3])
    print("Symbolic Jacobian:")
    sp.pprint(J)



# constants = dict()
# constants["Fin Moments"] = np.array([3,3,3])
# constants["correctiveConstants"] = 10
# constants["DragCoeff"] = 10
# constants["Mass"] = 10
# constants["Inertias"] = np.array([5,5,5])
# stateGuy = np.array([1,1,1,1,1,1])
# Jacobian = rocket.calculateANew(constants, stateGuy)

# print(Jacobian)





