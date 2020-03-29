import collections
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sklearn.linear_model import Ridge
import pandas as pd
from scipy.special import expit
import csv


class HillTypeMuscle:
    """
    Damped Hill-type muscle model adapted from Millard et al. (2013). The
    dynamic model is defined in terms of normalized length and velocity.
    To model a particular muscle, scale factors are needed for force, CE
    length, and SE length. These are given as constructor arguments.
    """

    def __init__(self):
        """
        :param Fmax: maximum isometric force
        :param k_shape: shape factor
        :param k_curve: curvature constant
        :param lam_ref: reference lambda
        :param lSE: length of SE element
        :param lSL: length of slack tendon
        :param lOpt: fascicle length
        """
        # INITIALIZE LATER CONSTANTS
        self.Fmax = 0 
        self.k_shape = 0
        self.k_curve = 0
        self.lam_ref = 0
        self.lSE = 0
        self.lSL = 0
        self.lOpt = 0

    def tendon_dynamics(self):
        # equations 7 & 8 in source [3]
        lam = (self.lSE - self.lSL)/self.lSL
        if lam <= 0:
            return 0
        else:
            numerator = np.exp((self.k_shape/self.lam_ref) * lam) - 1
            denominator = np.exp(self.k_shape) - 1
            return self.Fmax * (numerator/denominator)

    def get_length_CE(self):
        # returns CE dependent time
        return -1

    def get_velocity_CE(self):
        # returns velocity CE
        return -1

    def get_force_contractile_element(self, alpha, fL, fV):
        # source: [2], equation 1
        """
        :param alpha: FES activation (from our boi Trev)
        :param fL: active force length
        :param fV: active force velocity
        """
        return alpha*self.Fmax*fL*fV

    def get_force_length(self, w):
        return (-1/w**2)*((self.get_length_CE()/self.lOpt)**2) + (2/(w**2))*(self.get_length_CE()/self.lOpt) - 1/(w**2) + 1

    def get_force_velocity(self, f):
        #TODO: pass in Vmax
        Vmax = 0
        N = 1.5
        if self.get_velocity_CE() < 0:
            return (Vmax + self.get_velocity_CE()) / (Vmax - self.k_curve*self.get_velocity_CE())
        else:
            return N - ((N-1)*(Vmax + self.get_velocity_CE())) / (7.56*self.k_curve*self.get_velocity_CE() + Vmax)

    def get_force_parallel_element(self, w):
        if self.get_length_CE() <= (self.lOpt*(1-w)):
            return self.Fmax*((self.get_length_CE()-self.lOpt)/(self.lOpt*w))**2
        else:
            return 0

def plot_curves():
    """
    Plot force-length, force-velocity, SE, and PE curves.
    """
    lm = np.arange(0, 1.8, .01)
    vm = np.arange(-1.2, 1.2, .01)
    lt = np.arange(0, 1.07, .01)
    plt.subplot(2,1,1)
    plt.plot(lm, force_length_muscle(lm), 'r')
    plt.plot(lm, force_length_parallel(lm), 'g')
    plt.plot(lt, force_length_tendon(lt), 'b')
    plt.legend(('CE', 'PE', 'SE'))
    plt.xlabel('Normalized length')
    plt.ylabel('Force scale factor')
    plt.subplot(2, 1, 2)
    plt.plot(vm, force_velocity_muscle(vm), 'k')
    plt.xlabel('Normalized muscle velocity')
    plt.ylabel('Force scale factor')
    plt.tight_layout()
    plt.show()
