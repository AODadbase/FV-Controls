import numpy as np

class Kalman:
    def __init__(self, A, H, Q, R, P, x0):
        self.A = A
        self.H = H
        self.Q = Q
        self.R = R
        self.P = P
        self.x = x0

    def estimate(self):
        self.x = np.dot(self.A, self.x)
        self.P = np.dot(np.dot(self.A, self.P), self.A.T) + self.Q

    def update(self, z_k):
        inverted_part = np.dot(np.dot(self.H, self.P), self.H.T) + self.R
        K = np.dot(np.dot(self.P, self.H.T), np.linalg.inv(inverted_part))
        self.x = self.x + np.dot(K, z_k - np.dot(self.H, self.x))
        self.P = np.dot((np.eye(self.P.shape[0]) - np.dot(K, self.H)), self.P)
        
        


        
