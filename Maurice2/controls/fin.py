import numpy as np
from rocketpy.rocket.aero_surface.fins.trapezoidal_fins import TrapezoidalFins


class Fins(TrapezoidalFins):
    def __init__(
        self,
        n,
        root_chord,
        tip_chord,
        span,
        rocket_radius,
        cant_angle=0,
        sweep_length=None,
        sweep_angle=None,
        airfoil=None,
        name="Fins",
    ):
        super().__init__(
            n,
            root_chord,
            tip_chord,
            span,
            rocket_radius,
            cant_angle,
            sweep_length,
            sweep_angle,
            airfoil,
            name,
        )

        self.aileronAngles = np.array([0.0])

    def compute_forces_and_moments(
        self,
        stream_velocity,
        stream_speed,
        stream_mach,
        rho,
        cp,
        omega,
        *args,
    ): 
        R1, R2, R3, M1, M2, M3 = super().compute_forces_and_moments(
            stream_velocity,
            stream_speed,
            stream_mach,
            rho,
            cp,
            omega,
            *args,
        )

        M3 = M3 + self.computeAileronMoment(stream_velocity)
        otherM = self.computeOtherAileronMoment(stream_velocity)
        M2 = M2 + otherM[1]
        M1 = M1 + otherM[1]

        return R1, R2, R3, M1, M2, M3

    def computeAileronMoment(self, stream_velocity):
        delta1 = self.aileronAngles
        Mz = np.rad2deg(delta1[0])/8 * (-2.21e-09*(stream_velocity[2]**3) 
                            + 1.58e-06*(stream_velocity[2]**2) 
                            + 4.18e-06*stream_velocity[2] 
                            )
        return Mz

    def computeOtherAileronMoment(self):
        Mx = 0
        My = 0

        return [Mx, My]