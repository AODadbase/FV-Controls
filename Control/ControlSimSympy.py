from sympy import *
import numpy as np
from scipy import linalg
from control import lqr
import pandas as pd
import os

class Controls:
    def __init__(
            self,
            t_motor_burnout: float = 1.97,
        ):
        """Initialize the Controls class.

        Args:
            t_motor_burnout (float, optional): Time until motor burnout in seconds. Defaults to 1.97.
        """
        self.t_motor_burnout = t_motor_burnout # seconds
        self.csv_path = os.path.join(os.path.dirname(__file__), "important_data.csv")

    def getLineOfBestFitTime(self, motor_burnout: bool, var: str, n: int = 5):
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
        data = pd.read_csv(self.csv_path)
        t = data["# Time (s)"]
        y = None
        if (var == "mass"):
            y = data["Mass (g)"] / 1000  # Convert to kg
        elif (var == "inertia"):
            y = data["Longitudinal moment of inertia (kg·m²)"]
        elif (var == "drag force"):
            y = data["Drag force (N)"]
        elif (var == "drag coeff"):
            y = data["Drag coefficient (​)"]
        elif (var == "total vel"):
            y = data["Total velocity (m/s)"]
        elif (var == "angle of attack"):
            y = data["Angle of attack (°)"]
        else:
            raise ValueError(
                "Invalid variable. Choose from: " \
                "'mass'," \
                "'inertia', " \
                "'drag force', " \
                "'drag coeff', " \
                "'total vel', " \
                "'angle of attack'. "
                )

        # Filter data based on motor burnout
        if motor_burnout:
            mask = t >= self.t_motor_burnout
        else:
            mask = t < self.t_motor_burnout
        t = t[mask]
        y = y[mask]
        coeffs = np.polyfit(t, y, n)
        return coeffs, n

    def getLineOfBestFitAoA(self, var: str, n: int = 5):
        """Get the line of best fit for the given data with a polynomial of degree n.

        Args:
            var (str): The variable to fit the line to.
            n (int, optional): The degree of the polynomial to fit. Defaults to 5.

        Returns:
            tuple: A tuple containing the coefficients of the polynomial and its degree.
        """
        # Load the CSV data into a DataFrame
        data = pd.read_csv(self.csv_path)
        x = data["Angle of attack (°)"]
        if (var == "corrective moment coefficient"):
            y = data["Corrective moment Coefficient ()"]
        else:
            raise ValueError(
                "Invalid variable. Choose from: " \
                "'corrective moment coefficient'. "
                )

        mask = y != np.nan
        x = x[mask]
        y = y[mask]

        coeffs = np.polyfit(x, y, n)
        return coeffs, n

    def getConstants(self, t: float, AoA: float):
        """Get the constants for the rocket at time t.

        Args:
            t (float): The time in seconds.

        Returns:
            dict: A dictionary containing the fin moments, corrective constants, mass, inertia, drag coeff, and drag force.
        """

        constants = dict()
        I = Matrix([3, 3, 0.004]) # TODO: Check these values!!
        constants["fin moments"] = Matrix([0.003, 0.003, 0.003])  # Fin misalignment from open rocket

        motor_burnout = t > self.t_motor_burnout

        coeffs_mass, degree_mass = self.getLineOfBestFitTime(motor_burnout, "mass", self.csv_path)
        constants["mass"] = sum(coeffs_mass[i] * t**(degree_mass - i) for i in range(degree_mass + 1))
        
        coeffs_inertia, degree_inertia = self.getLineOfBestFitTime(motor_burnout, "inertia", self.csv_path)
        longI = sum(coeffs_inertia[i] * t**(degree_inertia - i) for i in range(degree_inertia + 1))
        I[0] = longI
        I[1] = longI
        constants["inertia"] = I
        
        coeffs_drag, degree_drag = self.getLineOfBestFitTime(motor_burnout, "drag coeff", self.csv_path)
        constants["drag coeff"] = sum(coeffs_drag[i] * t**(degree_drag - i) for i in range(degree_drag + 1))

        coeffs_drag_force, degree_drag_force = self.getLineOfBestFitTime(motor_burnout, "drag force", self.csv_path)
        constants["drag force"] = sum(coeffs_drag_force[i] * t**(degree_drag_force - i) for i in range(degree_drag_force + 1))


        # if t > self.t_motor_burnout:
        #     motor_burnout = True
        #     constants["mass"] = 2.589  # kg (empty motor)   

        #     I[0] = 0.345
        #     I[1] = 0.345

        #     coeffs_drag_coeff, degree_drag_coeff = self.getLineOfBestFitTime(motor_burnout, "drag coeff", self.csv_path)
        #     constants["drag coeff"] = sum(coeffs_drag_coeff[i] * t**(degree_drag_coeff - i) for i in range(degree_drag_coeff + 1))

        #     coeffs_drag_force, degree_drag_force = self.getLineOfBestFitTime(motor_burnout, "drag force", self.csv_path)
        #     constants["drag force"] = sum(coeffs_drag_force[i] * t**(degree_drag_force - i) for i in range(degree_drag_force + 1))

        # else:
        #     motor_burnout = False
        #     coeffs_mass, degree_mass = self.getLineOfBestFitTime(motor_burnout, "mass", self.csv_path)
        #     constants["mass"] = sum(coeffs_mass[i] * t**(degree_mass - i) for i in range(degree_mass + 1))
            
        #     coeffs_inertia, degree_inertia = self.getLineOfBestFitTime(motor_burnout, "inertia", self.csv_path)
        #     longI = sum(coeffs_inertia[i] * t**(degree_inertia - i) for i in range(degree_inertia + 1))
        #     I[0] = longI
        #     I[1] = longI
            
        #     coeffs_drag, degree_drag = self.getLineOfBestFitTime(motor_burnout, "drag coeff", self.csv_path)
        #     constants["drag coeff"] = sum(coeffs_drag[i] * t**(degree_drag - i) for i in range(degree_drag + 1))

        #     coeffs_drag_force, degree_drag_force = self.getLineOfBestFitTime(motor_burnout, "drag force", self.csv_path)
        #     constants["drag force"] = sum(coeffs_drag_force[i] * t**(degree_drag_force - i) for i in range(degree_drag_force + 1))

        #     coeffs_vel_mag, degree_vel_mag = self.getLineOfBestFitTime(motor_burnout, "total vel", self.csv_path)
        #     constants["total vel"] = sum(coeffs_vel_mag[i] * t**(degree_vel_mag - i) for i in range(degree_vel_mag + 1))

        # constants["inertia"] = I
        
        return constants

    def getAileronMoment(self, delta1: Symbol, v3: Symbol):
        """Get the aileron moment based on the aileron angle.

        Args:
            delta1 (float): The aileron angle in radians.

        Returns:
            Matrix: The symbolic moment vector [Mx, My, Mz] function of the aileron angle and rocket's vertical velocity.
        """
        
        M1, M2 = 0, 0
        # Verify later with data points
        M3 = delta1/8 * (4.11522634e-09*(v3**3) 
                        - 1.04938272e-06*(v3**2) 
                        + 3.35185185e-04*v3 
                        - 1.22222222e-02) # v3 = vertical velocity, Mz = roll moment

        return Matrix([M1, M2, M3])

    def getABC(self, t: float, AoA: float):
        """Get the equations of motion for the rocket.

        Args:
            t (float): The time in seconds.

        Returns:
            tuple: A tuple containing the A, B, and C Numpy arrays.
        """

        w1, w2, w3, v1, v2, v3 = symbols('w1 w2 w3 v1 v2 v3')
        I1, I2, I3 = symbols('I1 I2 I3')
        M1, M2, M3 = symbols('M1 M2 M3')
        mass, fg, rho, Cd, v_mag = symbols('m g rho Cd v_mag')
        delta1 = symbols('delta1')

        D = 1/2 * rho * Cd

        w1dot = ((I2 - I3) * w2 * w3 + M1) / I1
        w2dot = ((I3 - I1) * w3 * w1 + M2) / I2
        w3dot = ((I1 - I2) * w1 * w2 + M3) / I3
        v1dot = (D * v_mag * v1) / mass
        v2dot = (D * v_mag * v2) / mass
        v3dot = ((D * v_mag * v3) / mass) - fg

        f = Matrix([
            [w1dot],
            [w2dot],
            [w3dot],
            [v1dot],
            [v2dot],
            [v3dot]
        ])  
        g = Matrix([w3]) # Gyroscope measurement only measures w3 (w3 only relevant for roll stability)
        m = Matrix([w1, w2, w3, v1, v2, v3])
        n = Matrix([delta1])
        o = Matrix([w3])

        # Define parameters
        constants = self.getConstants(t, AoA)
        mass_rocket = constants["mass"]
        I = constants["inertia"]
        drag_coeff = constants["drag coeff"]
        Mf = constants["fin moments"]
        total_vel = constants["total vel"]
        Mf = self.getConstants(t)["fin moments"]
        M = Mf + self.getAileronMoment(delta1, v3)

        # C = D * v_mag / (sqrt(w1**2 + w2**2)) # Not included because might divide by 0

        params = {
            I1: Float(I[0]), # Ixx
            I2: Float(I[1]), # Iyy
            I3: Float(I[2]), # Izz
            M1: M[0],
            M2: M[1],
            M3: M[2],
            mass: Float(mass_rocket),
            fg: Float(9.81),
            rho: Float(1.225), # kg/m^3 temp constant rho
            Cd: Float(drag_coeff),
            v_mag: Float(total_vel),
        }

        delta1_e = 0.0  # equilibrium aileron angle, can be changed to optimize for different conditions
        f = f.subs(params)
        g = g.subs(params)
        m_e: dict = {
            w1: 0,
            w2: 0,
            w3: 0,
            v1: 0,
            v2: 0,
            v3: Float((mass * fg / (D * v_mag)).subs(params))  # equilibrium vertical velocity
        }
        n_e: dict = {delta1: delta1_e}

        A = f.jacobian(m).subs(m_e).subs(n_e)
        B = f.jacobian(n).subs(m_e).subs(n_e)
        C = g.jacobian(m).subs(m_e).subs(n_e)

        A = np.array(A).astype(np.float64)
        B = np.array(B).astype(np.float64)
        C = np.array(C).astype(np.float64)

        return A, B, C

def main():
    controls = Controls()
    A, B, C = controls.getABC(0)
    print("A Matrix:")
    print(A)
    print("B Matrix:")
    print(B)
    print("C Matrix:")
    print(C)

if __name__ == "__main__":
    main()

# # IGNORE TEMP CODE BELOW
# class ControlSimSympy:
#     def __init__(self, t: float):
#         """Initialize the ControlSimSympy class.

#         Args:
#             t (float): The time in seconds.
#         """
#         self.A, self.B, self.C = getABC(self, t)
#         self.Q = np.diag([1000, 1000, 1000, 10, 10, 10])  # State cost matrix
#         self.R = np.diag([1])  # Input cost matrix
#         self.K, self.S, self.E = lqr(self.A, self.B, self.Q, self.R)

#         self.A = np.array(self.A).astype(np.float64)
#         self.B = np.array(self.B).astype(np.float64)
#         self.C = np.array(self.C).astype(np.float64)
#         self.K = np.array(self.K).astype(np.float64)
#         self.Q = np.array(self.Q).astype(np.float64)
#         self.R = np.array(self.R).astype(np.float64)

#     def getControlInput(self, state: np.ndarray):
#         """Get the control input for the given state.

#         Args:
#             state (np.ndarray): The state vector.

#         Returns:
#             np.ndarray: The control input vector.
#         """
#         state = np.array(state).astype(np.float64)
#         control_input = -self.K @ state
#         return control_input