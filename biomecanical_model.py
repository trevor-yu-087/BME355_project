from scipy import interpolate
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# J \ddot{\theta} = F_{ta}d_{ta} - F_{gs} d_{ac} - F_{so}d_{ac} - F_g d_{cm}%0

def state_equation(t, x):
    """
    :param t: time
    :param x: state variable, [angle, angular velocity]
    :return: x_dot, [angular velocity, angular acceleration]
    """
    angular_velocity = x[1]
    angular_acceleration = (F_ta(t) * d_ta(x[0]) - F_gs(t) * d_at(x[0]) - F_so(t) * d_at(
        x[0]) - 1.05 * 9.81 * 0.05609 * np.cos(np.pi * x[0] / 180)) / 0.004248
    return [angular_velocity, angular_acceleration]


ta_time = np.arange(0, 20, .01)
ta_forces = np.sin(ta_time) * 40 + 40
gs_time = [0, 10]
gs_forces = [0, 10]
so_time = [0, 10]
so_forces = [0, 10]
TA_moment_arms = np.array([[-15.22613532, 0.05991416309],
                           [-0.1241071567, 0.04886266094],
                           [14.84963047, 0.0424248927],
                           [30.12968501, 0.04027896996]])
AT_moment_arms = np.array([[-15.02539273, 0.05385414642],
                           [0.1099468113, 0.06001116785],
                           [14.95690127, 0.06678542524],
                           [30.18731340, 0.07030750473]])

F_ta = interpolate.interp1d(ta_time, ta_forces)
F_gs = interpolate.interp1d(gs_time, gs_forces)
F_so = interpolate.interp1d(so_time, so_forces)
d_ta = interpolate.interp1d(TA_moment_arms[:, 0], TA_moment_arms[:, 1])
d_at = interpolate.interp1d(AT_moment_arms[:, 0], TA_moment_arms[:, 1])


def hit_max(t, x): return min(x[0] - 30, 0)
def hit_min(t, x): return max(x[0] + 15, 0)
hit_max.terminal = True
hit_min.terminal = True

t = 0
start = [0, 0]
time = []
solution = []
while(t < 10):
    sol = solve_ivp(state_equation, [t, 10], start, max_step=.001, events=(hit_max, hit_min))
    t = sol.t[-1] + .001
    start = [29.99 if sol.y[0, -1] > 0 else -14.99, 0]
    [time.append(t) for t in sol.t]
    [solution.append(x) for x in np.transpose(sol.y)]

solution = np.array(solution)
plt.subplot(3, 1, 1)
plt.plot(time, solution[:, 0])
plt.subplot(3, 1, 2)
plt.plot(time, solution[:, 1])
plt.subplot(3, 1, 3)
plt.plot(time, F_ta(time))
plt.plot(time, F_gs(time))
plt.show()
