import numpy as np
import matplotlib.pyplot as plt
from excitation_activation import FES_Activation, Fitted_Activation
import scipy
from scipy import integrate

class Hill_Type_Model:

    def __init__(self, muscle, alpha):
        """
        parameters chosen based on participant 3 in supplemntary materials (source [2])
        :param alpha: FES activation (from our boi Trev)
        :param Fmax: maximum isometric force
        :param k_shape: shape factor
        :param k_curve: curvature constant
        :param lam_ref: reference lambda
        :param lSE: length of SE element
        :param lOpt: fascicle length
        :param w: width of the active force length relation
        :param Vmax: maximum velocity
        :param p_angle: penation angle
        """
        self.alpha = alpha

        if(muscle == "Tibialis Anterior"):
            # CONSTANTS FOR TibAnt
            self.Fmax = 523
            self.k_shape = 2.73
            self.k_curve = 6.6
            self.lam_ref = 0.045
            self.lOpt = 0.21
            self.w = 0.49
            self.Vmax = 6
            self.p_angle = 5
        elif(muscle == "Soleus"):
            # CONSTANTS FOR Soleus
            self.Fmax = 4219
            self.k_shape = 2.51
            self.k_curve = 6.24
            self.lam_ref = 0.64
            self.lOpt = 0.266
            self.w = 0.8
            self.Vmax = 6.4
            self.p_angle = 25
        elif(muscle == "Gastrocnemius"):
            # CONSTANTS FOR gastrocnemius
            self.Fmax = 1816
            self.k_shape = 2.86
            self.k_curve = 8.4
            self.lam_ref = 0.051
            self.lOpt = 0.411
            self.w = 0.61
            self.Vmax = 4
            self.p_angle = 17

    def tendon_dynamics(self, lt):
        # equations 7 & 8 in source [3]
        lam = (lt - self.lOpt)/self.lOpt
        if lam <= 0:
            return 0
        else:
            numerator = np.exp((self.k_shape/self.lam_ref) * lam) - 1
            denominator = np.exp(self.k_shape) - 1
            return self.Fmax * (numerator/denominator)

    def get_force_length(self, lm):
        # source [2], equation 2
        return ((-1/self.w**2)*(lm/self.lOpt)**2 + (2/(self.w**2))*(lm/self.lOpt) - 1/(self.w**2) + 1)

    def get_force_velocity(self, vm):
        # source [2], equation 3
        N = 1.5
        if vm < 0:
            return (self.Vmax + vm)/(self.Vmax - self.k_curve*vm)
        else:
            return (N-((N-1)*(self.Vmax+vm)) / (7.56*self.k_curve*vm + self.Vmax))

    def get_force_contractile_element(self, lm, vm):
        # source: [2], equation 1
        return self.alpha*self.Fmax*lm*self.get_force_velocity(vm)

    def get_force_parallel_elastic(self,lt):
        # source [2], equation 4
        if lt <= (self.lOpt*(1-self.w)):
            return (self.Fmax*(lt-self.lOpt)/(self.lOpt*self.w))**2
        else: return 0

    def get_force_series_elastic(self, lm, lt, vm):
        # source: [3], equation adapted from simulation plan
        return ((self.get_force_contractile_element(lm,vm) + self.get_force_parallel_elastic(lt))*np.cos(self.get_angle(lm)))

    def get_angle(self, lm):
        # Source: [3] equation 11
        return np.arcsin((self.lOpt*np.sin(self.p_angle))/lm)

    def get_velocity_CE(self, lm, lt):
        sol = scipy.optimize.fsolve(self.model_dynamics, [0], (lm,lt))
        #print (sol)
        return sol

    def total_length_tendon(self, lm, lt):
        #Source: [3] equation 10
        return lt+lm*np.cos(self.get_angle(lm))


    def model_dynamics(self, vm, lm, lt):
        #Equation 9 source 3
        return self.get_force_series_elastic(lm,lt,vm)

    def simulate(self):

        def velocity_wrapper(t,x):
            return self.get_velocity_CE(x, self.total_length_tendon(x, self.lOpt))

        times = [0,2]
        solution = scipy.integrate.solve_ivp(velocity_wrapper, times, [5], rtol = 1e-8, atol = 1e-7)
        print(solution)
        plt.subplot(2,1,1)
        plt.plot(solution.t, solution.y.T)
        plt.ylabel('Normalized Length of CE')
        plt.subplot(2,1,2)
        plt.plot(solution.t, np.array(self.get_force_length(solution.y.T)))
        plt.ylabel('Force (N)')
        plt.xlabel('Time (s)')
        plt.show()

tibialis_anterior = Hill_Type_Model("Tibialis Anterior", 0.5)
tibialis_anterior.simulate()

gastrocnemius = Hill_Type_Model("Gastrocnemius", 0.2)
gastrocnemius.simulate()
