import numpy as np
import matplotlib.pyplot as plt
from excitation_activation import FES_Activation, Fitted_Activation
import scipy
from scipy import integrate

class Hill_Type_Model:

    def __init__(self, muscle, alpha, stim=None):
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


        # TODO figure out why these particular activations break it
        # self.alpha = lambda t: .2 - np.cos(t) * .2
        # self.alpha = lambda t: np.sin(t * 20) * .4 + .4

        # TODO figured it has something to do with an activation going from non-zero to zero.
        # TODO So I stopped alpha from becoming zero but hopefully will fix that eventually
        self.alpha = lambda t: max(alpha(t), .01)

        # all these work
        # self.alpha = lambda t: np.sin(t) * .2
        # self.alpha = lambda t: min(.5 * t, 1)
        # self.alpha = lambda t: 0 if t < 1 else 1

        if stim is None:
            self.stim = self.alpha
        else:
            self.stim = lambda t: max(stim(t), .01)

        if(muscle == "Tibialis Anterior"):
            # CONSTANTS FOR TibAnt
            self.Fmax = 523
            self.k_shape = 2.73
            self.k_curve = 6.6
            self.lam_ref = 0.045
            self.lOpt = 0.21  # m
            self.w = 0.49
            self.Vmax = 6
            self.p_angle = 5 * np.pi / 180
            self.lsl = 0.24
            self.FT = 25
            dm = 1.054  # g/cm3
            PCSA = 18.52    # cm2
            self.mass = 100*self.lOpt * dm * PCSA / 1000  # kg
        elif(muscle == "Soleus"):
            # CONSTANTS FOR Soleus
            self.Fmax = 4219
            self.k_shape = 2.51
            self.k_curve = 6.24
            self.lam_ref = 0.64
            self.lOpt = 0.266
            self.w = 0.8
            self.Vmax = 6.4
            self.p_angle = 25 * np.pi / 180
            self.lsl = 0.26
            self.FT = 20
            self.mass = None # not going to use this so we don't care that we don't have it
        elif(muscle == "Gastrocnemius"):
            # CONSTANTS FOR gastrocnemius
            self.Fmax = 1816
            self.k_shape = 2.86
            self.k_curve = 8.4
            self.lam_ref = 0.051
            self.lOpt = 0.411
            self.w = 0.61
            self.Vmax = 4
            self.p_angle = 17 * np.pi / 180
            self.lsl = .4
            self.FT = 50
            self.mass = None # this one is not used so we don't care that we don't have it

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
        if lm < 0: return 0
        return np.max(((-1/self.w**2)*(lm/self.lOpt)**2 + (2/(self.w**2))*(lm/self.lOpt) - 1/(self.w**2) + 1), 0)

    def get_force_velocity(self, vm):
        # source [2], equation 3
        N = 1.5
        if vm < 0:
            force = (self.Vmax + vm)/(self.Vmax - self.k_curve*vm)
        else:
            force =  (N-((N-1)*(self.Vmax+vm)) / (7.56*self.k_curve*vm + self.Vmax))
        return max(force, 0)

    def get_force_contractile_element(self, t, lm, vm):
        # source: [2], equation 1
        return max(self.alpha(t)*self.Fmax*self.get_force_length(lm)*self.get_force_velocity(vm), 0)

    def get_force_parallel_elastic(self,lm):
        # source [2], equation 4
        if lm <= (self.lOpt*(1-self.w)):
            return (self.Fmax*(lm-self.lOpt)/(self.lOpt*self.w))**2
        else: return 0

    def get_force_muscle(self, t, lm, vm):
        # source: [3], equation adapted from simulation plan
        return ((self.get_force_contractile_element(t, lm, vm) + self.get_force_parallel_elastic(lm))*np.cos(self.get_angle(lm)))

    def get_force_series_elastic_2(self, lm):
        # source: [3], equation 7
        result = self.Fmax * (np.exp(self.k_shape * (self.length_tendon(lm) - self.lsl) / (self.lam_ref * self.lsl)) - 1) / (np.exp(self.k_shape) - 1)
        return result

    def get_angle(self, lm):
        # Source: [3] equation 11
        return np.arcsin((self.lOpt*np.sin(self.p_angle))/lm)

    def get_velocity_CE(self, t, lm):
        sol = scipy.optimize.fsolve(self.model_dynamics, [0], (lm,t))
        #print("velocity")
        #print (sol)
        return sol

    def length_tendon(self, lm):
        #Source: [3] equation 10
        return self.lsl + self.lOpt - lm


    def model_dynamics(self, vm, lm, t):
        #Equation 9 source 3
        return self.get_force_series_elastic_2(lm) - self.get_force_muscle(t, lm, vm)

    def simulate(self, times, energy=False):

        def velocity_wrapper(t,x):
            vm = self.get_velocity_CE(t, x[0])
            if energy:
                x = [vm[0], self.E_dot(t, x[0], vm[0])]
                return x
            else:
                return vm

        if energy:
            initial_state = [self.lOpt + .001, 0]
        else:
            initial_state = [self.lOpt + .001]

        solution = scipy.integrate.solve_ivp(velocity_wrapper, times, initial_state, max_step=.001, rtol = 1e-8, atol = 1e-7)

        plt.subplot(4 if energy else 3,1,1)
        plt.plot(solution.t, solution.y.T[:, 0])
        plt.ylabel('Normalized Length of CE')
        plt.subplot(4 if energy else 3,1,2)
        plt.plot(solution.t, np.array(self.get_force_series_elastic_2(solution.y.T[:, 0])))
        plt.ylabel('Force (N)')
        plt.xlabel('Time (s)')
        plt.subplot(4 if energy else 3, 1, 3)
        plt.plot(solution.t, list(map(self.alpha, solution.t)))
        plt.ylabel("activation")
        if energy:
            plt.subplot(4, 1, 4)
            plt.plot(solution.t, solution.y.T[:, 1])
        plt.show()

    def E_dot(self, t, lm, vm):
        activation_maintenance = 0
        shorten_lengthen = 0
        mechanical_work = 0
        STIM = self.stim(t)
        ACT = self.alpha(t)
        h_AM = 1.28 * self.FT + 25
        A = STIM if STIM > ACT else (STIM + ACT) / 2
        S = 1.5
        F_ISO = self.get_force_length(lm) / lm
        if lm <= self.lOpt:
            activation_maintenance = h_AM * (A ** .6) * S
            if vm <= 0:
                shorten_lengthen = (-(vm/self.lOpt) * (100 - self.FT) / self.Vmax - 153 * (vm/self.lOpt) * self.FT / (100 * self.Vmax)) * (
                            A ** 2) * S
            else:
                shorten_lengthen = 400 * (vm/self.lOpt) * A * S / (vm/self.lOpt)
        else:
            activation_maintenance = (0.4 * h_AM + 0.6 * h_AM * F_ISO) * (A ** .6)
            if vm <= 0:
                shorten_lengthen = (-(vm/self.lOpt)* (100 - self.FT) / self.Vmax - 153 * (vm/self.lOpt) * self.FT / (
                            100 * self.Vmax)) * F_ISO * (A ** 2) * S
            else:
                shorten_lengthen = 400 * (vm/self.lOpt) * A * S * F_ISO / (vm/self.lOpt)

        mechanical_work = -(self.get_force_contractile_element(t, lm, vm) * vm) / self.mass
        return activation_maintenance + shorten_lengthen + mechanical_work


if __name__ == "__main__":
    tibialis_anterior = Hill_Type_Model("Tibialis Anterior", lambda t: np.sin(t * 20))
    tibialis_anterior.simulate([0, .5], energy=True)

    gastrocnemius = Hill_Type_Model("Gastrocnemius", lambda t: 0.2)
    gastrocnemius.simulate([0, .5])
