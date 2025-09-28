import numpy as np

from rocketpy.plots.aero_surface_plots import _TrapezoidalFinsPlots
from rocketpy.prints.aero_surface_prints import _TrapezoidalFinsPrints

from rocketpy.rocket.aero_surface.fins.trapezoidal_fins import TrapezoidalFins


class ourFins(TrapezoidalFins):
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

        self.aileronAngles = np.array([0,0,0,0])

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
        
        # print("The not roll moment is " + str(M2))
        #print("The roll moment is " + str(M3))

        return R1, R2, R3, M1, M2, M3

    def computeAileronMoment(self, stream_velocity):
        alphas = self.aileronAngles
        # 1000 * (alphas[1] + alphas[3] - alphas[2] - alphas [0])
        # return  0.02 * (alphas[1] + alphas[3] - alphas[2] - alphas [0]) * vScale
        Mz = alphas[0]/8 * (4.11522634e-09*(stream_velocity[2]**3) 
                            -1.04938272e-06*(stream_velocity[2]**2) 
                            + 3.35185185e-04*stream_velocity[2] 
                            -1.22222222e-02)
        return Mz


    def computeOtherAileronMoment(self, stream_velocity):
        vScale = np.sqrt(np.dot(stream_velocity,stream_velocity))/210

        alphas = self.aileronAngles
        # 1000 * (alphas[1] + alphas[3] - alphas[2] - alphas [0])
        Mx = 0
        My = 0
        
        return  [Mx, My]
