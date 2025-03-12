# def calculateA(self, time, velocities):
#         #Cols of A: 
#         #[w_1, w_2, w_3, v_x, v_y, v_z,]
        
#         index = self.times.index(time)
    
#         mass = self.masses[index]
        
#         Fg = mass * 9.81
#         thrust = self.thrusts[index]
#         drag_v = self.calculateDrag_v(time)
#         vMag = np.linalg.norm(velocities)

# #CHANGE!!!!!!!!!    
#         Cd = 10 #Get from open rocket
#         r = 1.225 #kg/m^3
#         Area = 10 #area of the diameter of the rocket body

#         VdotSclae = (-0.5 * Cd * r * Area * vMag / mass)
#         #theta dot = A theta --> make the eigenvalues negative for stability
#         #(F_d/m = v dot)
#         #F_d/m = deltaV

#         #corrective moment is the M term
#         #Have to divide moment by the right I

#         Cna = 10 #Normal Foce
#         dist = 10 #distance between center of mass and gravity
        
#         CMomentDir = vMag*vMag/((self.angular[0] * self.angular[0]  + self.angular[1] * self.angular[1]) **0.5 )
#         CMoment = CMomentDir * 0.5 * 1.225 * Area * (dist) * CMomentDir
#         #w1 : pitch
#         #w2 : yaw
#         #w3 : roll

#         A = [[0, -((self.inertias[2]-self.inertias[1])*self.angular[2])/self.inertias[0], -((self.inertias[2]-self.inertias[1])*self.angular[1])/self.inertias[0], 0, 0, 0],
#              [-((self.inertias[0]-self.inertias[2])*self.angular[2])/self.inertias[1], 0, -((self.inertias[0]-self.inertias[2])*self.angular[0])/self.inertias[1], 0, 0 ,0],
#              [-((self.inertias[1]-self.inertias[0])*self.angular[1])/self.inertias[2], -((self.inertias[1]-self.inertias[0])*self.angular[0])/self.inertias[2], 0, 0, 0, 0],
#              [0,0,0,VdotSclae,0,0],
#              [0,0,0,0,VdotSclae,0],
#              [0,0,0,0,0,VdotSclae - Fg/velocities[2]]]

#         return A
