# def E_dot():
#     activation_maintenance = 0
#     shorten_lengthen = 0
#     mechanical_work = 0
#     h_AM = 1.28 * FT + 25
#     A = STIM if STIM > ACT else (STIM + ACT) / 2
#     S = 1.5
#     F_ISO = get_force_length/get_length_CE
#     if L_CE <= L_CE_resting:
#         activation_maintenance = h_AM * (A**.6) * S
#         if V_CE_N <= 0:
#             shorten_lengthen = (-V_CE_N * (100 - FT)/V_CE_MAX_ST - 153 * V_CE_N * FT / (100 * V_CE_MAX_FT)) * (A**2) * S
#         else:
#             shorten_lengthen = 400 * V_CE_N * A * S / V_CE_N_ST
#     else:
#         activation_maintenance = (0.4 * h_AM + 0.6 * h_AM * F_ISO) * (A**.6)
#         if V_CE_N <= 0:
#             shorten_lengthen = (-V_CE_N * (100 - FT)/V_CE_MAX_ST - 153 * V_CE_N * FT / (100 * V_CE_MAX_FT)) * F_ISO * (A**2) * S
#         else:
#             shorten_lengthen = 400 * V_CE_N * A * S * F_ISO / V_CE_N_ST
#
#     mechanical_work = -(F_CE * V_CE)/m
#
#     return activation_maintenance * shorten_lengthen * mechanical_work