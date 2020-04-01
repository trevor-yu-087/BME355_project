from muscles import Hill_Type_Model
from excitation_activation import FES_Activation, Fitted_Activation
import matplotlib.pyplot as plt
import numpy as np
# from energy_expenditure import ...

if __name__ == '__main__':
    # todo: decide on a common time vector for all simulations, as some things need to be scaled
    time = 0

    """
    Make different activation profiles
    """
    ga_activation = Fitted_Activation('curve_datasets/gastrocnemius_activation.csv', width=0.06)
    # ga_activation.show_curves()
    sol_activation = Fitted_Activation('curve_datasets/soleus_activation.csv', width=0.09)
    # sol_activation.show_curves()

    """
    Make muscle objects for different simulations
    """
    # t = np.arange(0, 2, .01)
    # plt.plot(t, list(map(ga_activation.get_activation, t)))
    # plt.show()
    gastrocnemius = Hill_Type_Model("Gastrocnemius", ga_activation.get_activation)
    gastrocnemius.simulate([0, .8])

    """
    Simulate models with different activation profiles
    """

    """
    Run energy expenditure and physical dynamics models
    """
