from muscles import Hill_Type_Model
from excitation_activation import FES_Activation, Fitted_Activation
import matplotlib.pyplot as plt
import numpy as np
from regression import Regression
# from energy_expenditure import ...

def process_natural_stim(data_file):
    data = np.loadtxt(data_file, delimiter=',')
    norm_time = data[:, 0]/max(data[:,0])
    norm_amp = data[:, 1] / max(data[:, 1])
    # Extracted time from profile in which the envelope occurs
    t_profile = (0, 0.632)
    # Extracted start and end points of swing phase with respect to global walking cycle
    t_swing = (0.0456, 0.449)
    # Convert profile to time in global walking cycle
    wc_time = t_profile[1] * norm_time + 0.6
    wc_time = np.roll(wc_time, len(wc_time[wc_time > 1]))
    norm_amp = np.roll(norm_amp, len(wc_time[wc_time > 1]))
    wc_time[wc_time > 1] -= 1
    save = np.column_stack((wc_time, norm_amp))
    np.savetxt('curve_datasets/processed_natural_stimulation_2.csv', save, delimiter=',', fmt='%10.5f')


def get_fitted_natural_stimulation(time, scale=1):
    """
    :param scale: optional, how much to scale the voltage (normalized to max)
    :return: swing stance time, corresponding stimulation profile
    """
    data = np.loadtxt('curve_datasets/processed_natural_stimulation_2.csv', delimiter=',')
    centers = np.arange(0, 1, 0.005)
    model = Regression(data[:,0], data[:,1], centers, 0.04)
    # plt.plot(data[:,0], data[:,1])
    swing_stance_stim = model.eval(time)
    swing_stance_stim[swing_stance_stim < 0] = 0
    return scale*swing_stance_stim/max(swing_stance_stim)


def get_fitted_ankle_angle(time):
    data = np.loadtxt('curve_datasets/ankle_angle.csv', delimiter=',')
    centres = np.arange(0, 1, 0.005)
    sample_time = data[:, 0] / max(data[:, 0])
    sample_angle = data[:, 1] / max(data[:, 1])
    model = Regression(sample_time, sample_angle, centres, 0.04)
    swing_stance_aa = model.eval(time)
    return swing_stance_aa


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
    # ankle_angle = np.loadtxt('curve_datasets/ankle_angle.csv', delimiter=',')
    # aa_time = ankle_angle[:,0]/max(ankle_angle[:,0])
    # aa_angle = ankle_angle[:,1]/max(ankle_angle[:,1])
    # plt.plot(aa_time[aa_time > 0.6], aa_angle[aa_time > 0.6])
    angle = get_fitted_ankle_angle(time)
    plt.plot(time, angle, 'k')
    # For naturalistic stimulation, assume the maximum voltage is 41 V

    stim = get_fitted_natural_stimulation(time)
    plt.plot(time, stim)

    plt.legend(['gastroc act', 'sol act', 'regression angle', 'TA activation profile'])
    plt.show()

    # Two stimulation waveforms from paper
    # Constant rectangular
    fs = 40  # Hz
    pulse_width = 300  # micro seconds
    amplitude = 41  # V
    amp_rect = amplitude * np.ones(len(time))
    f_rect = fs * np.ones(len(time))
    rectangular_activation = FES_Activation(time, amp_rect, f_rect, 0.068, 0.076, 0.25)
    rectangular_activation.show_curves()

    # Enveloped naturalistic
    fs = 40  # Hz
    amplitude = 41  # V
    amp_enveloped = get_fitted_natural_stimulation(time, scale=amplitude)
    f_enveloped = fs * np.ones(len(time))
    enveloped_activation = FES_Activation(time, amp_enveloped, f_enveloped, 0.068, 0.076, 0.25)
    enveloped_activation.show_curves()

    """
    Make muscle objects for different simulations
    """
    # t = np.arange(0, 2, .01)
    # plt.plot(t, list(map(ga_activation.get_activation, t)))
    # plt.show()
    gastrocnemius = Hill_Type_Model("Gastrocnemius", ga_activation.get_activation)
    gastrocnemius.simulate([0, .8])
    tibialis = Hill_Type_Model("Tibialis Anterior", enveloped_activation.get_activation)
    tibialis.simulate([0.6, 1])
    """
    Simulate models with different activation profiles
    """

    """
    Run energy expenditure and physical dynamics models
    """
