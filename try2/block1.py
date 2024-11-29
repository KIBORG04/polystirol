from constants import *
from my_time import *
import math
import numpy as np
import plotly.graph_objects as go
import scipy.integrate as sci

def block1(n):
    return kimptAndKfc(n)

class Rolf:
    y = np.array(())
    t = np.array(())

def RKSolver(f, time, y0, args):
    t_values = [time[0]]
    y_values = [y0]
    t = time[0]
    y = np.array(y0)
    h = dt
    while t < time[1]:
        if t + h > time[1]:  # Корректируем последний шаг, чтобы не превысить t_end
            h = time[1] - t

        k1 = h * np.array(f(t, y, *args))
        k2 = h * np.array(f(t + h/2, y + k1/2, *args))
        k3 = h * np.array(f(t + h/2, y + k2/2, *args))
        k4 = h * np.array(f(t + h, y + k3, *args))

        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
        t = t + h

        t_values.append(t)
        y_values.append(y.tolist())

    rol = Rolf()
    rol.y = np.array(y_values).T
    rol.t = time

    return rol

def kimptAndKfc(n):
    y0 = [0, 0, T0[n]] + [0 for _ in range(FACTIONS)]
    sol = sci.solve_ivp(kimpt_pend, [timeStart, timeEnd], y0, args=(n,), t_eval=time)
    #sol = RKSolver(kimpt_pend, [timeStart, timeEnd], y0, args=#(n,))
    return sol

def kimpt_pend(t, y, n):
    Xi, Xs, T = y[0], y[1], y[2]
    P = y[3:]

    # Константы кинетики
    Ki = Ki0 * math.pow(math.e, -Ei/(R*T))
    Ks = Ks0 * math.pow(math.e, -Es/(R*T))
    Cp = (Ms[n] * Xs) / Mr[n]
    Kf = A - B * math.pow(math.e, C * Cp)

    dXi = Ki * ms0[n] * (1 - Xi) * (1 - Xs)
    dXs = (Ki * Jk[n] * (1 - Xi) + Ks * math.pow(Jk[n] * Xi, 0.5 + ws*Cp))*(1-Xs)

    numeratorLeft = (Ktes * dXs * ms0[n] * Vt[n])
    numeratorRigth = (Kf * Fst * L[n] * (T - Th[n]) * Ghs[n] * Ch * ph)/(Kf * L[n] * Fst + Ghs[n] * Ch * ph)
    denominator = (Map * Cap) + (Vt[n] * dpm * Cpm)
    dT = (numeratorLeft - numeratorRigth) / denominator

    # Константы концентрации
    ks = Ks / (math.pow(Jk[n], 0.5-ws*Cp))
    ki = Ki
    J = Jk[n] * (1 - Xi)
    if Xi == 0:
        J = 0
    ms = ms0[n] * (1 - Xi)
    dP = [0 for _ in range(FACTIONS)]

    for i in range(FACTIONS):
        if i == 0:
            dp = -ks * ms * P[i] + ki * ms * J
        elif i == FACTIONS-1:
            dp = ks * ms * P[i-1]
        else:
            dp = -ks * ms * P[i] + ks * ms * P[i-1]
        dP[i] = dp

    dydt = [dXi,
            dXs,
            dT] + dP
    return dydt

def main():
    block1(0)

if __name__ == "__main__":
    main()
