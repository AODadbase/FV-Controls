from sympy import *
from sympy.algebras.quaternion import Quaternion
import numpy as np
from scipy import linalg
from control import lqr
import pandas as pd
import os

class Controls:
    def __init__(
            self,
            t_motor_burnout: float = 1.971,
        ):
        """Initialize the Controls class. Rocket body axis is aligned with y-axis.

        Args:
            t_motor_burnout (float, optional): Time until motor burnout in seconds. Defaults to 1.971.
        """
        self.t_motor_burnout = t_motor_burnout # seconds
        self.csv_path = os.path.join(os.path.dirname(__file__), "important_data.csv")

    def getLineOfBestFitTime(self, var: str, n: int = 1):
        """Get the line of best fit for the given data with a polynomial of degree n.

        Args:
            var (str): The variable to fit the line to.
            n (int, optional): The degree of the polynomial to fit. Defaults to 1.

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
        else:
            raise ValueError(
                "Invalid variable. Choose from: " \
                "'mass'," \
                "'inertia'. " \
                )

        # Filter data based on motor burnout
        mask = t >= self.t_motor_burnout
        t = t[mask]
        y = y[mask]
        coeffs = np.polyfit(t, y, n)
        return coeffs, n


    def getLineOfBestFitAoA(self, var: str, n: int = 5):
        """Get the line of best fit for the given data with a polynomial of degree n. Choose from "stability margin" or "normal force coeff".

        Args:
            var (str): The variable to fit the line to.
            n (int, optional): The degree of the polynomial to fit. Defaults to 5.

        Returns:
            tuple: A tuple containing the coefficients of the polynomial and its degree.
        """
        # Load the CSV data into a DataFrame
        data = pd.read_csv(self.csv_path)
        x = data["Angle of attack (°)"]
        if (var == "stability margin"):
            y = data["Stability margin calibers (​)"]
        elif (var == "normal force coeff"):
            y = data["Normal force coefficient (​)"]
        else:
            raise ValueError(
                "Invalid variable. Choose from: " \
                "'stability margin', " \
                "'normal force coeff'. "
                )

        mask = y != np.nan
        motor_burnout = data["# Time (s)"] > self.t_motor_burnout # TODO: use motor_burnout to split data between pre and post burnout
        x = x[mask]
        y = y[mask]

        coeffs = np.polyfit(x, y, n)
        return coeffs, n

    def getLineOfBestFitVel(self, var: str, n: int = 2):
        """Get the line of best fit for the given data with a polynomial of degree n.

        Args:
            var (str): The variable to fit the line to.
            n (int, optional): The degree of the polynomial to fit. Defaults to 2.

        Returns:
            tuple: A tuple containing the coefficients of the polynomial and its degree.
        """
        # Load the CSV data into a DataFrame
        # TODO: set data to __init__ for efficiency
        data = pd.read_csv(self.csv_path)
        x = data["Total velocity (m/s)"]
        if (var == "drag force"):
            y = data["Drag force (N)"]
        elif (var == "drag coeff"):
            y = data["Drag coefficient (​)"]
        else:
            raise ValueError(
                "Invalid variable. Choose from: " \
                "'drag force', " \
                "'drag coeff'. "
                )
        motor_burnout = data["# Time (s)"] > self.t_motor_burnout 
        x = x[not motor_burnout]
        y = y[not motor_burnout]

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
        I = Matrix([0.004, 0.004, 0.29]) # Post burnout inertia values from OpenRocket, kg*m^2
        m = 2.589  # Post burnout mass from OpenRocket, kg
        T = Matrix([0., 0., 0.])  # N

        motor_burnout = t > self.t_motor_burnout

        # TODO: for added efficiency, only call getLineOfBestFitTime once per variable and store the results
        if not motor_burnout:
            coeffs_mass, degree_mass = self.getLineOfBestFitTime("mass")
            m = sum(coeffs_mass[i] * t**(degree_mass - i) for i in range(degree_mass + 1))

            coeffs_inertia, degree_inertia = self.getLineOfBestFitTime("inertia")
            I[2] = sum(coeffs_inertia[i] * t**(degree_inertia - i) for i in range(degree_inertia + 1))

            times = pd.read_csv(self.csv_path)["# Time (s)"]
            thrust = pd.read_csv(self.csv_path)["Thrust (N)"]
            T[2] = np.interp(t, times, thrust) # Thrust acting in z direction

        constants["inertia"] = I
        constants["mass"] = m
        constants["thrust"] = T
        constants["fin moments"] = Matrix([0.003, 0.003, 0.003])  # Fin misalignment from open rocket
        
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

    def R_BW_from_q(self, qw, qx, qy, qz):
        """Convert a quaternion to a rotation matrix. World to body frame.

        Args:
            qw (float): The scalar component of the quaternion.
            qx (float): The x component of the quaternion.
            qy (float): The y component of the quaternion.
            qz (float): The z component of the quaternion.

        Returns:
            Matrix: The rotation matrix from world to body frame.
        """
        s = (qw**2 + qx**2 + qy**2 + qz**2)**-Rational(1,2) # Normalizing factor
        qw, qx, qy, qz = qw*s, qx*s, qy*s, qz*s # Normalized quaternion components

        xx,yy,zz = qx*qx, qy*qy, qz*qz
        wx,wy,wz = qw*qx, qw*qy, qw*qz
        xy,xz,yz = qx*qy, qx*qz, qy*qz
        return Matrix([
            [1-2*(yy+zz),   2*(xy+wz),   2*(xz-wy)],
            [2*(xy-wz),     1-2*(xx+zz), 2*(yz+wx)],
            [2*(xz+wy),     2*(yz-wx),   1-2*(xx+yy)]
        ])

    def getAB(self, t: float, AoA: float):
        """Get the equations of motion for the rocket.

        Args:
            t (float): The time in seconds.
            AoA (float): The angle of attack in radians.

        Returns:
            tuple: A tuple containing the A and B Numpy arrays.
        """

        w1, w2, w3, v1, v2, v3 = symbols('w1 w2 w3 v1 v2 v3', real = True) # Angular and linear velocities
        qw, qx, qy, qz = symbols('qw qx qy qz', real = True) # Quaternion components
        I1, I2, I3 = symbols('I1 I2 I3', real = True) # Moments of inertia
        M1, M2, M3 = symbols('M1 M2 M3', real = True) # Moments
        T1, T2, T3 = symbols('T1 T2 T3', real = True) # Thrusts
        mass, rho, A, g = symbols('m rho A g', real = True) # Mass, air density, reference area, gravity
        delta1 = symbols('delta1', real = True) # Aileron angle

        eps = Float(1e-10)  # Small term to avoid division by zero
        v = Matrix([v1, v2, v3]) # Velocity vector
        v_mag = v.norm() + eps # Magnitude of velocity with small term to avoid division by zero
        vhat = v / v_mag  # Unit vector in direction of velocity

        ## Thrust ##
        T : Matrix = Matrix([T1, T2, T3])  # Thrust vector, T1 and T2 are assumed 0

        ## Gravity ##
        Fg_world = Matrix([0.0, 0.0, -mass * g])
        R_world_to_body = self.R_BW_from_q(qw, qx, qy, qz)  # Rotation matrix from world to body frame
        Fg : Matrix = R_world_to_body * Fg_world  # Transform gravitational force to body frame

        ## Drag Force ##
        Fd_mag = -(0.627 + -0.029*v_mag + 1.95e-3*v_mag**2) # Drag force approximation
        Fd : Matrix = Fd_mag * vhat # Drag force vector

        ## Lift Force ##
        beta = Piecewise((atan2(v2, v1), v1**2 + v2**2 > eps),
                 (Float(0), True))
        L = 1/2 * rho * v_mag**2 * (2 * pi * AoA) * A # Lift force approximation
        nL = Matrix([
            -cos(AoA) * cos(beta),
            -cos(AoA) * sin(beta),
            sin(AoA)
        ]) # Lift direction unit vector
        Fl : Matrix = L * nL # Lift force vector

        ## Total Forces ##
        F = T + Fd + Fl + Fg # Total force vector

        S = Matrix([[0, -w3, w2],
                    [w3, 0, -w1],
                    [-w2, w1, 0]])
        q_vec = Matrix([qw, qx, qy, qz])
        Omega = Matrix([
            [0, -w1, -w2, -w3],
            [w1, 0, w3, -w2],
            [w2, -w3, 0, w1],
            [w3, w2, -w1, 0]
        ])

        ## Corrective Moment Coefficient ##

        Cn = ... # TODO: implement this with a sinusoidal fit, copy paste function here like for SM
        Cnalpha = Cn / AoA
        SM = 2.8 + -0.48*AoA + 0.163*AoA**2 + -0.0386*AoA**3 + 5.46E-03*AoA**4 + -4.61E-04*AoA**5 + 2.28E-05*AoA**6 + -6.1E-07*AoA**7 + 6.79E-09*AoA**8
        d = 7.87/100 # m, rocket diameter
        C = v_mag**2 * A * Cnalpha * (SM * d) * rho / 2

        w1dot = ((I2 - I3) * w2 * w3 - C * w1 + M1) / I1
        w2dot = ((I3 - I1) * w3 * w1 - C * w2 + M2) / I2
        w3dot = ((I1 - I2) * w1 * w2 + M3) / I3
        vdot = F/mass - S * v
        qdot = (Omega * q_vec) * Float(1/2)

        f = Matrix([
            [w1dot],
            [w2dot],
            [w3dot],
            [vdot[0]],
            [vdot[1]],
            [vdot[2]],
            [qdot[0]],
            [qdot[1]],
            [qdot[2]],
            [qdot[3]]
        ])

        m = Matrix([w1, w2, w3, v1, v2, v3, qw, qx, qy, qz])
        n = Matrix([delta1])

        # Define parameters
        constants = self.getConstants(t, AoA)
        mass_rocket = constants["mass"]
        I = constants["inertia"]
        T = constants["thrust"]
        Mf = constants["fin moments"]
        M = Mf + self.getAileronMoment(delta1, v3)

        params = {
            I1: Float(I[0]), # Ixx
            I2: Float(I[1]), # Iyy
            I3: Float(I[2]), # Izz
            M1: M[0],
            M2: M[1],
            M3: M[2],
            T1: T[0],
            T2: T[1],
            T3: T[2],
            mass: Float(mass_rocket),
            rho: Float(1.225), # kg/m^3 temp constant rho
            A: Float(pi*(7.87/100/2)**2), # m^2 reference area
            g: Float(9.81), # m/s^2
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

        A = np.array(A).astype(np.float64)
        B = np.array(B).astype(np.float64)

        return A, B
    
    def getC(self, t: float, AoA: float):
        """Get the C matrix for the rocket at time t.

        Args:
            t (float): The time in seconds.
            AoA (float): The angle of attack in radians.

        Returns:
            np.ndarray: The C matrix.
        """
        
# main() doesn't work yet, just a placeholder for now for future testing
def main():
    controls = Controls()
    A, B = controls.getAB(0, 0)
    print("A Matrix:")
    print(A)
    print("B Matrix:")
    print(B)

if __name__ == "__main__":
    main()
