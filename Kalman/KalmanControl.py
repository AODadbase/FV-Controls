import numpy as np
from Equations import Kalman
from Equations import TrueValues
from eulerEquations import state
import matplotlib.pyplot as plt
import random


class Controller:
    def __init__(self, initialPosition, initialVelocity):
        self.pos = initialPosition
        self.velo = initialVelocity
    