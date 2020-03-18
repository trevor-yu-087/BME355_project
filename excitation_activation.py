import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt


class FES_Activation:
    """
    Activation-Excitation dynamics of an artificially activated muscle from Romero et al. (2015)
    f_stim is a stimulation frequency profile
    u_stim is a stimulation intensity profile
    t is the time vector associated with the profiles
    """
    def __init__(self, time, u_stim, f_stim, t_rise, t_fall, FT_percent):
        self.u_stim = u_stim    #need to set one of these to a vector of 1s
        self.f_stim = f_stim
        self.t_rise = t_rise
        self.t_fall = t_fall
        self.time = time
        # Params used by Romero et al. from other literature sources
        self.a2 = 2.5
        self.R = 40
        self.rCF = (self.R - 91.2) / (-1.03)
        self.f0 = self.R*np.log((self.a2-1)*np.exp(self.rCF/self.R)-self.a2)
        self.a1 = -self.a2*np.exp(-self.f0/self.R)
        self.u_threshold = (5+((70-1)*(50-5))/((127-1)))     # measured by Romero et al.
        self.u_saturation = (5+((110-1)*(50-5))/((127-1)))   # measured by Romero et al.
        self.set_Te(FT_percent)

    def set_Te(self, FT_percent):
        """
        Sets excitation time constant based on fast twitch fibres composition and other constants from literature
        :param FT_percent: percentage fast twitch fibres in the muscle
        """
        FT = 2549 #
        S0 = 450000 # [N / cm ^ 2]
        gamma = 10 * np.pi / 180 # [rad]
        dm = 1054 # [kg / m ^ 3]
        lF = 0.1 # [m]

        PCSA = FT / (S0 * np.cos(gamma))
        m = PCSA * lF * dm
        self.Te = (25 + 0.1 * m * (1 - FT_percent)) / 1000

    def Sf(self,f_stim):
        """
        :param f_stim: frequency stimulation profile
        :return: frequency scaling factor for each point in the profile
        """
        return (self.a1 - self.a2)/(1 + np.exp((f_stim - self.f0)/self.R)) + self.a2

    def Su(self, u_stim):
        """
        Intensity scaling can have either pulse width or amplitude as profile
        :param u_stim: intensity stimulation profile
        :return: intensity scaling factor for each point in the profile
        """
        # Su as function of time
        Su = []
        for u in u_stim:
            if u < self.u_threshold:
                Su.append(0)
            elif self.u_threshold <= u < self.u_saturation:
                Su.append((u - self.u_threshold)/(self.u_saturation - self.u_threshold))
            else:
                Su.append(1)
        return np.array(Su)

    def get_excitation_signal(self, f_stim, u_stim):
        return self.Sf(f_stim) * self.Su(u_stim)

    def excitation_finite_difference(self, t, ex):
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
        forward = np.interp(forward_time, self.time, ex)
        backward = np.interp(backward_time, self.time, ex)
        # forward = self.elastance(forward_time)
        # backward = self.elastance(backward_time)
        return (forward - backward) / (2*dt)

    def get_activation_signal(self, time, f_stim, u_stim):
        sol = integrate.solve_ivp(self.activation_derivative, (0, max(time)), [0, 0], max_step=0.01)
        return sol.t, sol.y


    def activation_derivative(self, t, a):
        """
        Activation derivative for the ODE of the 2nd block in Hammerstein model
        :param t: time
        :param a: [a, a_dot]
        :return: [a_dot, a_dot_dot]
        """
        ex = self.get_excitation_signal(self.f_stim, self.u_stim)
        if self.excitation_finite_difference(t, ex) > 0:
            k1 = self.Te * self.t_rise
            k2 = self.Te + self.t_rise
        else:
            k1 = self.Te * self.t_fall
            k2 = self.Te + self.t_fall
        A = [[0,1],[-1/k1, -k2/k1]]
        ex_int = np.interp(t, self.time, ex)
        return np.matmul(A, a) + ex_int*np.array([0, 1/k1])


if __name__ == "__main__":
    time = np.linspace(0, 2, 100)
    f_stim = np.ones(100)
    u_stim = np.zeros(100)
    u_stim[25:49] = 2
    u_stim[75:99] = 2
    t_rise = 0.068 #[s]
    t_fall = 0.076 #[s]
    TA_Activation = FES_Activation(time, u_stim, f_stim, t_rise, t_fall, 0.25)
    t,y = TA_Activation.get_activation_signal(time, f_stim, u_stim)
    plt.figure()
    plt.plot(t,y[0,:])
    plt.show()
