import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from regression import Regression


class FES_Activation:
    """
    Activation-Excitation dynamics of an artificially activated muscle from Romero et al. (2015)
    f_stim is a stimulation frequency profile. Centre frequency is 39.6 Hz.
    u_stim is a stimulation intensity profile. Note the threshold intensity is 29 and saturation of 43
    time is the time vector associated with the profiles
    t_rise and t_fall are time constants for the muscle activation
    FT_percent is the percentage of fast-twitch fibres in the muscle
    I might add params for the cross sectional area
    """
    def __init__(self, time, u_stim, f_stim, t_rise, t_fall, FT_percent):
        self.u_stim = u_stim
        self.f_stim = f_stim
        self.t_rise = t_rise  # sec
        self.t_fall = t_fall  # sec
        self.FT_percent = FT_percent
        self.time = time
        # Params used by Romero et al. from other literature sources
        self.a2 = 2.5
        self.R = 40  # 15 ?
        self.rCF = (self.R - 91.2) / (-1.03)
        self.f0 = self.R*np.log((self.a2-1)*np.exp(self.rCF/self.R)-self.a2)
        self.a1 = -self.a2*np.exp(-self.f0/self.R)
        self.u_threshold = (5+((70-1)*(50-5))/((127-1)))     # measured by Romero et al.
        self.u_saturation = (5+((110-1)*(50-5))/((127-1)))   # measured by Romero et al.
        self.Te = 0.025     # Estimate for entire muscle excitation time
        self.ex = self.get_excitation_signal()
        t, y = self.get_activation_signal()
        self.t_act = t  # For interpolation
        self.y_act = y[0, :]    # For interpolation

    def Sf(self):
        """
        :return: frequency scaling factor for each point in the frequency profile
        """
        freq_scale = (self.a1 - self.a2)/(1 + np.exp((self.f_stim - self.f0)/self.R)) + self.a2
        freq_scale[freq_scale > 1] = 1
        return freq_scale

    def Su(self):
        """
        Intensity scaling can have either pulse width or amplitude as profile
        :param u_stim: intensity stimulation profile
        :return: intensity scaling factor for each point in the profile
        """
        Su = []
        for u in self.u_stim:
            if u < self.u_threshold:
                Su.append(0)
            elif self.u_threshold <= u < self.u_saturation:
                Su.append((u - self.u_threshold)/(self.u_saturation - self.u_threshold))
            else:
                Su.append(1)
        return np.array(Su)

    def get_excitation_signal(self):
        return self.Sf() * self.Su()

    def excitation_finite_difference(self, t):
        """
        Calculates finite-difference approximation of excitation derivative for determining whether to use
        rise or fall time constants.
        :param t: time
        :param ex: excitation signal
        :return: finite-difference approximation of time derivative of excitation signal
        """
        dt = .01
        forward_time = t + dt
        backward_time = max(0, t - dt) # small negative times are wrapped to end of cycle
        forward = np.interp(forward_time, self.time, self.ex)
        backward = np.interp(backward_time, self.time, self.ex)
        return (forward - backward) / (2*dt)

    def get_activation_signal(self):
        sol = integrate.solve_ivp(self.activation_derivative, (min(self.time), max(self.time)), [0, 0], max_step=0.01)
        return sol.t, sol.y

    def get_activation(self, t):
        """
        :param t: time to evaluate activation signal at
        :return: linear interpolation of activation signal given by solve_ivp
        """
        return np.interp(t, self.t_act, self.y_act)

    def get_excitation(self, t):
        """
        :param t: time to evaluate excitation signal at
        :return: linear interpolation of excitation signal
        """
        return np.interp(t, self.time, self.ex)

    def activation_derivative(self, t, a):
        """
        Activation derivative for the ODE of the 2nd block in Hammerstein model
        :param t: time
        :param a: [a, a_dot]
        :return: [a_dot, a_dot_dot]
        """
        if self.excitation_finite_difference(t) > 0:
            k1 = self.Te * self.t_rise
            k2 = self.Te + self.t_rise
        else:
            k1 = self.Te * self.t_fall
            k2 = self.Te + self.t_fall
        A = [[0, 1], [-1/k1, -k2/k1]]
        ex_int = np.interp(t, self.time, self.ex)
        return np.matmul(A, a) + ex_int*np.array([0, 1/k1])

    def show_curves(self):
        t, y = self.get_activation_signal()
        plt.figure()
        plt.plot(t,y[0,:])
        # plt.plot(self.time, self.get_activation(self.time))
        plt.plot(self.time, self.u_stim/self.u_saturation)
        plt.plot(self.time, self.get_excitation_signal())
        plt.legend(['Activation', 'FES (normalized)', 'Excitation'])
        plt.show()


class Fitted_Activation:
    def __init__(self, data_file, width=0.04):
        """
        Constructor requires activation fit data
        :param data_file: path to .csv file containing activation data
        :param width: optional, set Gaussian width
        """
        data = np.loadtxt(data_file, delimiter=',')
        self.walking_cycle_percent = data[:, 0]
        self.activation = data[:, 1]
        self.time = np.arange(0, 1, 0.002)
        self.centers = np.arange(0, 1, 0.005)
        self.width = width
        self.model = Regression(self.walking_cycle_percent, self.activation, self.centers, self.width)

    def show_curves(self):
        plt.figure()
        plt.plot(self.walking_cycle_percent, self.activation, 'k')
        plt.plot(self.time, self.get_activation(self.time), 'r')
        plt.xlabel('Walking Cycle (%)')
        plt.ylabel('Activation')
        plt.legend(['Data', 'Regression'])
        plt.show()

    def get_activation(self, t):
        """
        :param t: time point at which to evaluate activation, in percentage of walking cycle
        :return: activation level at that time
        """
        a = self.model.eval(t)
        a[a < 0.002] = 0
        return a


if __name__ == "__main__":
    # Example FES Activation
    time = np.linspace(0, 2, 100)
    f_stim = 66*np.ones(100)
    u_stim = np.zeros(100)
    u_stim[25:49] = 40
    u_stim[75:99] = 40
    # U between 29 and 43
    # F0 = 39.6 Hz
    t_rise = 0.068  # [s]
    t_fall = 0.076  # [s]
    FT_percent = 0.25
    TA_Activation = FES_Activation(time, u_stim, f_stim, t_rise, t_fall, FT_percent)
    # TA_Activation.show_curves()

    # Using fitted activation
    ga_activation = Fitted_Activation('curve_datasets/gastrocnemius_activation.csv', width=0.06)
    # ga_activation.show_curves()
    # print(ga_activation.get_activation(0.5))

    sol_activation = Fitted_Activation('curve_datasets/soleus_activation.csv', width=0.09)
    # sol_activation.show_curves()

    data = np.loadtxt('curve_datasets/artificial_stimulation_regimen.csv', delimiter=',')
    plt.plot(data[:,0], data[:,1])
    plt.show()


