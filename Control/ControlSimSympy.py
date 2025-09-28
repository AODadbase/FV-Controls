from sympy import *
import numpy as np
from scipy import linalg
from control import lqr
import pandas as pd
import os

t_motor_burnout = 1.97  # seconds

def getLineOfBestFit(motor_burnout, var, csv_path, n=5):
    """Get the line of best fit for the given data with a polynomial of degree n.

    Args:
        motor_burnout (bool): If True, use data after motor burnout. If False, use data before motor burnout.
        var (str): The variable to fit the line to.
        csv_path (str): The path to the CSV file containing historical data.
        n (int, optional): The degree of the polynomial to fit. Defaults to 5.

    Returns:
        tuple: A tuple containing the coefficients of the polynomial and its degree.
    """
    # Load the CSV data into a DataFrame
    data = pd.read_csv(csv_path)
    t = data["# Time (s)"]
    y = None
    if (var == "mass"):
        y = data["Mass (kg)"]
    elif (var == "inertia"):
        y = data["Longitudinal moment of inertia (kg·m²)"]
    elif (var == "drag force"):
        y = data["Drag force (N)"]
    elif (var == "drag coeff"):
        y = data["Drag coefficient"]
    else:
        raise ValueError("Invalid variable. Choose from 'mass', 'inertia', 'drag force', or 'drag coeff'.")

    # Filter data based on motor burnout
    if motor_burnout:
        mask = t >= t_motor_burnout
    else:
        mask = t < t_motor_burnout
    t = t[mask]
    y = y[mask]
    coeffs = np.polyfit(t, y, n)
    return coeffs, n


def getConstants(self, t):
    """Get the constants for the rocket at time t.

    Args:
        t (float): The time in seconds.

    Returns:
        dict: A dictionary containing the fin moments, corrective constants, mass, inertia, drag coeff, and drag force.
    """
    # Define the path to the CSV file
    csv_path = os.path.join(os.path.dirname(__file__), "../Dynamics/MAURICE2_openrocket.csv")

    constants = dict()
    I = np.array([3, 3, 0.004])
    constants["fin moments"] = np.array([0.003, 0.003, 0.003])  # Fin misalignment from open rocket
    constants["corrective constants"] = 1

    if t > t_motor_burnout:
        constants["mass"] = 2.589  # kg (empty motor)

        I[0] = 0.345
        I[1] = 0.345

        coeffs, degree = getLineOfBestFit(True, "drag coeff", csv_path)
        constants["drag coeff"] = sum(coeffs[i] * t**(degree - i) for i in range(degree + 1))

        coeffs_drag_force, degree_drag_force = getLineOfBestFit(True, "drag force", csv_path)
        constants["drag force"] = sum(coeffs_drag_force[i] * t**(degree_drag_force - i) for i in range(degree_drag_force + 1))

    else:
        coeffs_mass, degree_mass = getLineOfBestFit(False, "mass", csv_path)
        constants["mass"] = sum(coeffs_mass[i] * t**(degree_mass - i) for i in range(degree_mass + 1))
        
        coeffs_inertia, degree_inertia = getLineOfBestFit(False, "inertia", csv_path)
        longI = sum(coeffs_inertia[i] * t**(degree_inertia - i) for i in range(degree_inertia + 1))
        I[0] = longI
        I[1] = longI
        
        coeffs_drag, degree_drag = getLineOfBestFit(False, "drag coeff", csv_path)
        constants["drag coeff"] = sum(coeffs_drag[i] * t**(degree_drag - i) for i in range(degree_drag + 1))

        coeffs_drag_force, degree_drag_force = getLineOfBestFit(False, "drag force", csv_path)
        constants["drag force"] = sum(coeffs_drag_force[i] * t**(degree_drag_force - i) for i in range(degree_drag_force + 1))

    constants["inertia"] = I
    
    return constants

def EoM(self, time):
    """Get the equations of motion for the rocket.

    Args:
        time (float): The time in seconds.
    """
    constants = getConstants(time)
    m = constants["mass"]
    I = constants["inertia"]
    Cd = constants["drag coeff"]
    D = constants["drag force"]
    F = constants["fin moments"]
    Mk = constants["corrective constants"]

    w1, w2, w3, v1, v2, v3 = symbols('w1 w2 w3 v1 v2 v3')

    C = D * (v1**2 + v2**2 + v3**2) / (sqrt(w1**2 + w2**2))

    x = Matrix([w1, w2, w3, v1, v2, v3])

