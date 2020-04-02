from scipy import interpolate
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from muscles import Hill_Type_Model


# J \ddot{\theta} = F_{ta}d_{ta} - F_{gs} d_{ac} - F_{so}d_{ac} - F_g d_{cm}%0

TA_moment_arms = np.array([[-15.22613532, 0.05991416309],
                           [-0.1241071567, 0.04886266094],
                           [14.84963047, 0.0424248927],
                           [30.12968501, 0.04027896996]])
AT_moment_arms = np.array([[-15.02539273, 0.05385414642],
                           [0.1099468113, 0.06001116785],
                           [14.95690127, 0.06678542524],
                           [30.18731340, 0.07030750473]])

d_ta = interpolate.interp1d(TA_moment_arms[:, 0], TA_moment_arms[:, 1])
d_at = interpolate.interp1d(AT_moment_arms[:, 0], TA_moment_arms[:, 1])

def simulate(F_ta, F_gs, F_so, start, times, plot=False):
    def state_equation(t, x):
        """
        :param t: time
        :param x: state variable, [angle, angular velocity]
        :return: x_dot, [angular velocity, angular acceleration]
        """
        print("start")
        print(t)
        print(x)
        angular_velocity = x[1]
        angular_acceleration = (F_ta(t) * d_ta(x[0]) - F_gs(t) * d_at(x[0]) - F_so(t) * d_at(
            x[0]) - 1.05 * 9.81 * 0.05609 * np.cos(np.pi * x[0] / 180)) / 0.004248
        return [angular_velocity, angular_acceleration]

    def hit_max(t, x): return min(x[0] - 29, 0)
    def hit_min(t, x): return max(x[0] + 14, 0)
    hit_max.terminal = True
    hit_min.terminal = True

    time = []
    solution = []
    while(times[0] < times[1]):
        sol = solve_ivp(state_equation, times, start, max_step=.001, events=(hit_max, hit_min))
        times[0] = sol.t[-1] + .001
        start = [29 if sol.y[0, -1] > 0 else -14, 0]
        [time.append(t) for t in sol.t]
        [solution.append(x) for x in np.transpose(sol.y)]

    solution = np.array(solution)
    if plot:
        plt.subplot(3, 1, 1)
        plt.plot(time, solution[:, 0])
        plt.ylabel("ankle angle degrees")
        plt.subplot(3, 1, 2)
        plt.plot(time, solution[:, 1])
        plt.ylabel("ankle angular velocity degrees/second")
        plt.subplot(3, 1, 3)
        plt.plot(time, F_ta(time))
        plt.plot(time, F_gs(time))
        plt.plot(time, F_so(time))
        plt.ylabel("force N")
        plt.xlabel("time seconds")
        plt.legend(["Tibialis Anterior", "Gastrocnemius", "Soleus"])
        plt.show()

    return solution


if __name__ == "__main__":
    ta = Hill_Type_Model("Tibialis Anterior", lambda t: t)
    gs = Hill_Type_Model("Gastrocnemius", lambda t: 1 - t)

    ta_sol, ta_force = ta.simulate([0, 1], plot=True)
    gs_sol, gs_force = gs.simulate([0, 1], plot=True)

    print(len(ta_force))
    print(len(gs_force))

    F_ta = interpolate.interp1d(ta_sol.t, ta_force)
    F_gs = interpolate.interp1d(gs_sol.t, gs_force)
    F_so = interpolate.interp1d([0, 1], [500, 200])

    simulate(F_ta, F_gs, F_so, [0, 0], [0, 1], plot=True)
