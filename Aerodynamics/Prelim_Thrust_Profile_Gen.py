class Prelim_Thrust_Profile_Gen:
    def __init__(self, m_dot, p_exhaust_velocity, p_exhaust, p_atm, area_exhaust_exit):
    self.m_dot = 2 # kg/s
    self.exhaust_velocity = 2500 #m/s
    p_exhaust = .5 #atm, exhaust pressure
    p_atm = 3 # lolz
    area_exhaust_exit = .01 #m^2
    def thrust_gen():
        Force_Thrust = m_dot * exhaust_velocity + (p_exhaust - p_atm)*area_exhaust_exit