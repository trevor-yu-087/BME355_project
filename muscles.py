import numpy as np
import matplotlib.pyplot as plt
import scipy

class Hill_Type_Model:

    def __init__(self, muscle, alpha, resting_length_muscle, resting_length_tendon):
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
        self.resting_length_muscle = resting_length_muscle
        self.resting_length_tendon = resting_length_tendon

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

    def norm_tendon_length(self, muscle_tendon_length, normalized_muscle_length):
        """
        :param muscle_tendon_length: non-normalized length of the full muscle-tendon
            complex (typically found from joint angles and musculoskeletal geometry)
        :param normalized_muscle_length: normalized length of the contractile element
            (the state variable of the muscle model)
        :return: normalized length of the tendon
        """
        return (muscle_tendon_length - self.resting_length_muscle * normalized_muscle_length) / self.resting_length_tendon

    def tendon_dynamics(self, lt):
        # equations 7 & 8 in source [3]
        lam = (self.norm_tendon_length(lt, 1) - self.lOpt)/self.lOpt
        if lam <= 0:
            return 0
        else:
            numerator = np.exp((self.k_shape/self.lam_ref) * lam) - 1
            denominator = np.exp(self.k_shape) - 1
            return self.Fmax * (numerator/denominator)

    def get_force_length(self, lt):
        # source [2], equation 2
        return (-1/self.w**2)*((self.norm_tendon_length(lt, 1)/self.lOpt)**2) + (2/(self.w**2))*(self.norm_tendon_length(lt,1)/self.lOpt) - 1/(self.w**2) + 1

    def get_force_velocity(self, lm, lt):
        # source [2], equation 3
        vm = []
        N = 1.5
        for i in range(len(lm)):
            if self.get_velocity_CE(lm, lt) < 0:
                vm.append(self.Vmax + self.get_velocity_CE(lm, lt)[i]) / (self.Vmax - self.k_curve*self.get_velocity_CE(lm, lt)[i])
            else:
                vm.append(N - ((N-1)*(self.Vmax + self.get_velocity_CE(lm, lt))) / (7.56*self.k_curve*self.get_velocity_CE(lm, lt) + self.Vmax))
        return np.asarray(vm)

    def get_force_contractile_element(self, lm, lt):
        # source: [2], equation 1
        fCE = []
        for i in range(len(lm)):
            fCE.append(self.alpha*self.Fmax*self.get_force_length(lt[i])*self.get_force_velocity(lm[i],lt[i]))
        return np.asarray(fCE)

    def get_force_parallel_elastic(self,lm):
        # source [2], equation 4
        fPE = []
        for i in range(len(lm)):
            if self.norm_tendon_length(lm[i], 1) <= (self.lOpt*(1-self.w)):
                fPE.append(self.Fmax*((self.norm_tendon_length(lm[i], 1)-self.lOpt)/(self.lOpt*self.w))**2)
            else:
                fPE.append(0)
        return np.asarray(fPE)

    def get_force_series_elastic(self, lm, lt):
        # source: [3], equation adapted from simulation plan
        fSE = []
        for i in range(len(lm)):
            fSE.append((self.get_force_contractile_element(lm[i],lt[i]) + self.get_force_parallel_elastic(lm[i]))*np.cos(self.p_angle))
        return np.asarray(fSE)

    def get_velocity_CE(self, lm, lt):
        """
        :param a: activation (between 0 and 1)
        :param lm: normalized length of muscle (contractile element)
        :param lt: normalized length of tendon (series elastic element)
        :return: normalized lengthening velocity of muscle (contractile element)
        """
        beta = 0.1 # damping coefficient (see damped model in Millard et al.)
        sol = scipy.optimize.fsolve(self.model_dynamics, [0], (self.alpha,lm,lt,self.Fmax,beta))
        #print (sol)
        return sol

    def model_dynamics(self, vm, lm, lt):
        return self.alpha*self.get_force_length(lt)*self.get_force_velocity(lm,lt) + self.get_force_parallel_elastic(lm) + 0.1*vm-self.get_force_length(lt)

    def simulate(self):
        length = self.resting_length_tendon+self.resting_length_muscle
        
        def velocity_wrapper(t,x):
            return self.get_velocity_CE(x, self.norm_tendon_length(length, x))

        ICs = [0,2]
        solution = scipy.integrate.solve_ivp(velocity_wrapper, ICs, [1], rtol = 1e-8, atol = 1e-7)

        normal_tendon_length = self.norm_tendon_length(length, solution.y[0])
        plt.subplot(2,1,1)
        plt.plot(solution.t, solution.y.T)
        plt.ylabel('Normalized Length of CE')
        plt.subplot(2,1,2)
        plt.plot(solution.t, np.array(self.get_force_length(normal_tendon_length))*self.Fmax)
        plt.ylabel('Force (N)')
        plt.xlabel('Time (s)')
        plt.show()

def plot_curves():
    """
    Plot force-length, force-velocity, SE, and PE curves.
    """
    # lm = np.arange(0, 1.8, .01)
    # vm = np.arange(-1.2, 1.2, .01)
    # lt = np.arange(0, 1.07, .01)
    # plt.subplot(2,1,1)
    # plt.plot(lm, force_length_muscle(lm), 'r')
    # plt.plot(lm, force_length_parallel(lm), 'g')
    # plt.plot(lt, force_length_tendon(lt), 'b')
    # plt.legend(('CE', 'PE', 'SE'))
    # plt.xlabel('Normalized length')
    # plt.ylabel('Force scale factor')
    # plt.subplot(2, 1, 2)
    # plt.plot(vm, force_velocity_muscle(vm), 'k')
    # plt.xlabel('Normalized muscle velocity')
    # plt.ylabel('Force scale factor')
    # plt.tight_layout()
    # plt.show()
