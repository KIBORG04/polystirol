from constants import *
from my_time import *
import math
import scipy.integrate as sci

def block2(n, P):
    return kimptAndKfc(n, P)

def kimptAndKfc(n, P):
    y0 = [0, T1[n]] + P
    sol = sci.solve_ivp(kimpt_pend, [timeStart, timeEnd], y0, args=(n,), t_eval=time)
    return sol

def kimpt_pend(t, y, n):
    Xd, T = y[0], y[1]
    P = y[2:]

    # Константы кинетики
    cps = Ms[n] / (Md1[n] + Ms[n] + Mr[n])
    Cp = cps + (Md1[n] * Xd) / Mr[n]
    Kf = A - B * math.pow(math.e, C * Cp)
    Kd = Kd0 * (1 - b * Cp) * math.pow(math.e, -Ed / (R * T))

    dXd = (Kd * math.pow(Jacd1[n], 0.25 + wd*Cp)) * (1 - Xd)

    numeratorLeft = (Kted * dXd * md10[n] * Vt[n])
    numeratorRigth = (Kf * Fst * L[n] * (T - Th[n]) * Ghd1[n] * Ch * ph)/(Kf * L[n] * Fst + Ghd1[n] * Ch * ph)
    denominator = (Map * Cap) + (Vt[n] * dpm * Cpm)
    dT = (numeratorLeft - numeratorRigth) / denominator

    # Константа концентрации
    kd = Kd/math.pow(Jacd1[n], 0.75-wd*Cp)
    md = md10[n]*(1-Xd)
    dP = [0 for _ in range(FACTIONS)]
    for i in range(FACTIONS):
        if i == 0:
            dp = -kd * md * P[i]
        elif i == FACTIONS-1:
            dp = kd * md * P[i-1]
        else:
            dp = -kd * md * P[i] + kd * md * P[i-1]
        dP[i] = dp

    dydt = [
        dXd,
        dT] + dP
    return dydt
