from muscles import Hill_Type_Model
from excitation_activation import FES_Activation, Fitted_Activation
import matplotlib.pyplot as plt
import numpy as np
from regression import Regression
# from energy_expenditure import ...

def process_natural_stim(data_file):
    data = np.loadtxt(data_file, delimiter=',')
    norm_time = data[:, 0]/max(data[:,0])
    scaled_amp = data[:, 1] / max(data[:, 1])
    # Extracted time from profile in which the envelope occurs
    t_profile = (0, 0.632)
    # Extracted start and end points of swing phase with respect to global walking cycle
    t_swing = (0.0456, 0.449)
    # Convert profile to time in global walking cycle
    wc_time = t_profile[1] * norm_time
    save = np.column_stack((wc_time, scaled_amp))
    np.savetxt('curve_datasets/processed_natural_stimulation.csv', save, delimiter=',', fmt='%10.5f')
    # Fit regression

def get_fitted_natural_stimulation(scale=1):
    """
    :param scale: optional, how much to scale the voltage (normalized to max)
    :return: swing stance time, corresponding stimulation profile
    """
    data = np.loadtxt('curve_datasets/processed_natural_stimulation.csv', delimiter=',')
    centers = np.arange(0, 0.55, 0.005)
    model = Regression(data[:,0], data[:,1], centers, 0.04)
    # plt.plot(data[:,0], data[:,1])
    # Select amplitude values that occur during swing
    swing_stance_time = np.linspace(0, 0.4, 200)
    swing_stance_stim = model.eval(swing_stance_time)
    return 0.6 + swing_stance_time, scale*swing_stance_stim


if __name__ == '__main__':
    duration = 0.4  # seconds, based on 40% of walking cycle at 1 Hz

    """
    Make different activation profiles
    """
    ga_activation = Fitted_Activation('curve_datasets/gastrocnemius_activation.csv', width=0.06)
    # ga_activation.show_curves()
    sol_activation = Fitted_Activation('curve_datasets/soleus_activation.csv', width=0.09)
    # sol_activation.show_curves()

    time = np.linspace(0.6, 1, 300)
    plt.plot(time, ga_activation.get_activation(time))
    plt.plot(time, sol_activation.get_activation(time))
    ankle_angle = np.loadtxt('curve_datasets/ankle_angle.csv', delimiter=',')
    aa_time = ankle_angle[:,0]/max(ankle_angle[:,0])
    aa_angle = ankle_angle[:,1]/max(ankle_angle[:,1])
    plt.plot(aa_time[aa_time > 0.6], aa_angle[aa_time > 0.6])

    # For naturalistic stimulation, assume the maximum voltage is 41 V

    t, stim = get_fitted_natural_stimulation()
    plt.plot(t, stim)

    plt.legend(['gastroc act', 'sol act', 'angle', 'TA activation profile'])
    plt.show()


    """
    Make muscle objects for different simulations
    """

    """
    Simulate models with different activation profiles
    """

    """
    Run energy expenditure and physical dynamics models
    """
