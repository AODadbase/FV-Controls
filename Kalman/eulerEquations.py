import numpy as np
import csv
from scipy.optimize import approx_fprime
import sympy as sp

class state:
    def __init__(self, angular, inertias, torques, filePath):
        #Define all the class variables
        self.angular = np.array(angular)
        self.inertias = np.array(inertias)
        self.torques = np.array(torques)
        self.filePath = filePath

    def parseData(self):
        #take data from rocket and extract angular velocity, inertial components, and torque components
        self.times = []
        self.verticalVelocities = []
        self.masses = []
        self.thrusts = []
        self.dragForces = []
        self.dragCoefficients = []
        self.airTemperatures = []
        self.airPressures = []
        with open(self.filePath, "r", newline="") as file:
            reader = csv.reader(file, delimiter=",")
            for row in reader:
                if row[0][0] != '#':
                    self.times.append(row[0])
                    self.verticalVelocities.append(row[1])
                    self.masses.append(row[2])
                    self.thrusts.append(row[3])
                    self.dragForces.append(row[4])
                    self.dragCoefficients.append(row[5])
                    self.airTemperatures.append(row[6])
                    self.airPressures.append(row[7])
            return self.times, self.verticalVelocities, self.masses, self.thrusts, self.dragForces, self.dragCoefficients, self.airTemperatures, self.airPressures
    
    def calculateDrag_v(self, time):
        GASCONSTANT = 287.052874
        index = self.times.index(time)
        verticalVelocity = self.verticalVelocities[index]
        dragForce = self.dragForces[index]
        dragCoefficient = self.dragCoefficients[index]
        airTemperature = self.airTemperatures[index]
        airPressure = self.airPressures[index]
        
        airDensity = airPressure/(airTemperature*GASCONSTANT)
        area = (2*dragForce)/(airDensity*dragCoefficient*(verticalVelocity**2))
        drag_v = airDensity*verticalVelocity*dragCoefficient*area
        return drag_v


    def aLambdaFunc(self, x, FinMoments,Mk, Dk , mass ,inertias):
        #x = (w1, w2, w3, v1, v2, v3)
        velo = np.array([x[3], x[4], x[5]])
        wPrime = np.array([x[0], x[1]])


        correctiveCoeff = Mk * np.dot(velo, velo) / np.linalg.norm(wPrime)
        vmag = np.linalg.norm(velo)

        w1Dot = (1/inertias[0]) *((inertias[2] - inertias[1]) * x[2]*x[1] + (correctiveCoeff * x[0]) + FinMoments[0])
        w2Dot = (1/inertias[1]) * ((inertias[2] - inertias[0])*x[0]*x[2] + (correctiveCoeff*x[1]) + FinMoments[1])
        w3Dot = (1/inertias[2]) * ((inertias[0] - inertias[1]) *x[0]*x[1] + FinMoments[2])
        v1Dot = (Dk * vmag * x[3])/mass
        v2Dot = (Dk * vmag * x[4])/mass
        v3Dot = (Dk * vmag * x[5])/mass - 9.8

        return np.array([w1Dot, w2Dot,w3Dot, v1Dot, v2Dot, v3Dot ])
    
    def calculateANew(self, constants, state):
        FinMoments = constants["Fin Moments"]
        Mk = constants["correctiveConstants"]
        Dk = constants["DragCoeff"]
        mass = constants["Mass"]
        inertias = constants["Inertias"]
        
        physicsWrapped = lambda x : self.aLambdaFunc(x,FinMoments, Mk, Dk, mass, inertias)
        epsilon = 1e-6
        J = approx_fprime(state, physicsWrapped, epsilon)
        return J

    
    def calculateMoments(velocities, alphas):
        #Fill in once amruth gives me STUFF
        return np.array([1,2,3])

    def calculateBFunctional(x, alphas, inertias, self):
        moments = self.calculateMoments(np.array([x[3], x[4], x[5]]), alphas)
        w1Dot = moments[0]/inertias[0]
        w2Dot = moments[1]/inertias[1]
        w3Dot = moments[2]/inertias[2]

        return (w1Dot, w2Dot, w3Dot, 0,0,0)
    
    def calculateBNew(alphas, constants, self):
        inputsWrapped = lambda x : self.calculateBFunctional(x, alphas, constants["Inertias"], self)
        epsilon = 1e-6
        defaultAngles = np.array([0,0,0,0])
        J = approx_fprime(defaultAngles, inputsWrapped, epsilon)
        return J

    def getState(self, time):
        index = self.times.index(time)
        verticalVelocity = self.verticalVelocities[index]
        
        currentState = [
            self.angular[0],
            self.angular[1],
            self.angular[2],
            0,
            0,
            verticalVelocity
            ]
        return currentState

        

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



def jacobian():
    angular = [0, 0, 0]
    torques = [0, 0, 0]
    inertias = [0, 0, 0]
    filePath = "flightVehicleDataV3.csv"

    rocket = state(angular, torques, inertias, filePath)
    symPyJacobian()

# constants = dict()
# constants["Fin Moments"] = np.array([3,3,3])
# constants["correctiveConstants"] = 10
# constants["DragCoeff"] = 10
# constants["Mass"] = 10
# constants["Inertias"] = np.array([5,5,5])
# stateGuy = np.array([1,1,1,1,1,1])
# Jacobian = rocket.calculateANew(constants, stateGuy)

# print(Jacobian)






