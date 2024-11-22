import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math


# Масса аппарата (реактора), кг
Map = 16840
# Теплоёмкость хладагента, кДж/(кг К)
Ch = 4.19
# Плотность хладагента, кг/л
Ph = 1
# Теплоёмкость материала аппарата, Дж/(кгК)
Cap = 0.46
# Молекулярная масса стирола, кг/моль
MMs = 0.092
# Молекулярная масса дивинила, кг/моль
MMd = 0.054
# Площадь поверхности теплообмена в реакторе, m^2
Fst = 42.0
# Плотность реакционной массы, кг/л
dpm = 0.763
# Теплоёмкость реакционной массы, кДж/(кгК)
Cpm = 1.8436
# энергия активации инициирования Ei [Дж/моль ]
Ei = 59460
# энергия активации роста полистирольных цепей Es [Дж/моль ]
Es = 71184
#
Ed = 80850
# универсальная газовая постоянная, [кДж/(К моль)]
R = 8.32
# постоянные констант скорости инициирования цепи [л/(мин моль)]
Ki0 = 0.835 * math.pow(10, 10)
#
k = 2 * math.pow(10, 12)
#
b = 3.77
# постоянные констант скорости роста цепи
Ks0 = 5.76 * math.pow(10, 11)
# шаг дискретизация
dt = 0.01
# коэффициенты идентификации [Дж/(м2*К*мин)]
A = 14.17
# Коэффициент тепловедение стирола
Kted = 80.9
# Коэффициент тепловедение дивинила
Ktes = 74.87
# Длительность процесса
tn = 30
B = 1.007
C = 10.38

# параметр зависимости, учитывающей влияние концентрации
# полимера на порядок реакции по катализатору
ws = -0.4142
wd = -0.184
eps = 0.0001

# количество фракций полимера
n1 = 3

# Начальные температуры
T0 = np.array([310.4, 309.6, 312.7, 309.3, 313.0, 311.2])
# Концентрация катализатора Jk , моль/л.
Jk = np.array([0.00478, 0.0045, 0.00416, 0.00361, 0.0041, 0.00402])
# Концентрация катализатора по данным дозировки Jkdoz,моль/л
Jkdoz = np.array([0.0983, 0.00955, 0.00894, 0.00822, 0.0079, 0.00919])
# Концентрация стирола ms0, моль/л.
ms0 = np.array([0.745, 0.7304, 0.7129, 0.7106, 0.6683, 0.7168])
# Температура хладагента Тh, °С
Th = np.array([290, 289.5, 288, 290, 296, 288.5])
# Расход хладагента Gh, л/мин
Gh = np.array([1500, 2000, 2000, 1666, 1500, 1500])
Ms = np.array([1140, 1080, 1085, 1081, 998, 1085])

# Динамическая вязкость
ng = np.array([0.125, 0.129, 0.13, 0.147, 0.129, 0.137])
Mr = np.array([11400, 11140, 11370, 11500, 11270, 11290])
Md1 = np.array([1765, 1687, 1701, 1664, 1558, 1702])
Md2 = np.array([837, 715, 715, 722, 663, 710])
pr = np.array([0.748, 0.748, 0.739, 0.75, 0.744, 0.739])
ps = np.array([0.91, 0.916, 0.92, 0.919, 0.921, 0.92])
pd = np.array([0.638, 0.637, 0.637, 0.635, 0.637, 0.637])
Td0 = np.array([307.4, 310.4, 309.5, 313.3, 314.5, 309.4])

Vk = 139.2

#
Vt = np.array([(Vk + ((Ms[i] / ps[i]) + (Mr[i] / pr[i]))) for i in range(len(Ms))])


# расчет степени заполнения
L = np.array([((Vt[i] - 3000) / 28500) for i in range(len(Vt))])


#
Vd = np.array([(Md1[i] / pd[i]) for i in range(len(pd))])


#
Jac = np.array([(Jk[i] * ((Vt[i]) / (Vt[i] + Vd[i]))) for i in range(len(Vt))])


#
md0 = np.array([((pd[i] * Vd[i]) / (MMd * (Vt[i] + Vd[i]))) for i in range(len(pd))])


def kimpt1(Xi, Xs, T, n):
    cp = (Ms[n] * Xs) / Mr[n]
    Kf = A - (B * math.e ** (C * cp))
    Ki = Ki0 * math.pow(math.e, -(Ei / (R * T)))
    Ks = Ks0 * (math.e ** -(Es / (R * T)))

    Xi1 = Ki * (ms0[n]) * (1 - Xi) * (1 - Xs)

    Xs1 = (Ki * Jk[n] * (1 - Xi)) + (
        Ks * ((Jk[n] * Xi) ** (0.5 + (ws * cp))) * (1 - Xs)
    )

    T1 = (
        (Ktes * Xs1 * ms0[n] * Vt[n])
        - (
            (Kf * Fst * L[n] * (T - Th[n]) * Gh[n] * Ch * Ph)
            / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
        )
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    Xi2 = Ki * (ms0[n]) * (1 - (Xi + Xi1 * dt / 2)) * (1 - (Xs + Xs1 * dt / 2))

    Xs2 = (Ki * Jk[n] * (1 - Xi + Xi1 * dt / 2)) + (
        Ks
        * ((Jk[n] * (Xi + Xi1 * dt / 2)) ** (0.5 + (ws * cp)))
        * (1 - Xs + Xs1 * dt / 2)
    )

    T2 = (
        (Ktes * (Xs2) * ms0[n] * Vt[n])
        - (Kf * Fst * L[n] * ((T + T1 * dt / 2) - Th[n]) * Gh[n] * Ch * Ph)
        / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    Xi3 = Ki * (ms0[n]) * (1 - (Xi + Xi2 * dt / 2)) * (1 - (Xs + Xs2 * dt / 2))
    Xs3 = (Ki * Jk[n] * (1 - Xi + Xi2 * dt / 2)) + (
        Ks
        * ((Jk[n] * (Xi + Xi2 * dt / 2)) ** (0.5 + (ws * cp)))
        * (1 - Xs + Xs2 * dt / 2)
    )

    T3 = (
        (Ktes * Xs3 * ms0[n] * Vt[n])
        - (Kf * Fst * L[n] * ((T + T2 * dt / 2) - Th[n]) * Gh[n] * Ch * Ph)
        / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    Xi4 = Ki * (ms0[n]) * (1 - (Xi + Xi3 * dt)) * (1 - (Xs + Xs3 * dt))
    Xs4 = (Ki * Jk[n] * (1 - Xi + Xi3 * dt / 2)) + (
        Ks
        * ((Jk[n] * (Xi + Xi2 * dt / 2)) ** (0.5 + (ws * cp)))
        * (1 - Xs + Xs3 * dt / 2)
    )

    T4 = (
        (Ktes * Xs4 * ms0[n] * Vt[n])
        - (Kf * Fst * L[n] * ((T + T3 * dt) - Th[n]) * Gh[n] * Ch * Ph)
        / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    xi = Xi + (dt * (Xi1 + 2 * (Xi2 + Xi3) + Xi4) / 6)

    xs = Xs + (dt * (Xs1 + 2 * (Xs2 + Xs3) + Xs4) / 6)
    t = T + dt * (T1 + 2 * (T2 + T3) + T4) / 6
    ms = ms0[n] * (1 - xs)
    J = Jk[n] * (1 - xi)
    ks = Ks / (math.pow(Jk[n], (0.5 - (ws * cp))))
    ki = Ki

    return xi, xs, t, ms, J, ks, ki


def kfc1(ms, J, ks, ki, p):
    pp = np.zeros(n1)
    mms = np.zeros(n1)
    P1 = 0
    P2 = 0
    P3 = 0
    P4 = 0
    for i in range(n1):
        if i == 0:
            P1 = ki * ms * J - ks * ms * p[i]
            P2 = ki * ms * J - ks * ms * (p[i] + P1 * dt / 2)
            P3 = ki * ms * J - ks * ms * (p[i] + P2 * dt / 2)
            P4 = ki * ms * J - ks * ms * (p[i] + P3 * dt / 2)

        elif (i > 0) & (i < n1 - 2):
            P1 = (-ks * ms * p[i]) + (ks * ms * p[i - 1])
            P2 = (-ks * ms * (p[i] + P1 * dt / 2)) + (ks * ms * p[i - 1] + P1 * dt / 2)
            P3 = (-ks * ms * (p[i] + P2 * dt / 2)) + (ks * ms * p[i - 1] + P2 * dt / 2)
            P4 = (-ks * ms * (p[i] + P3 * dt / 2)) + (ks * ms * p[i - 1] + P3 * dt / 2)

        else:

            P1 = ks * ms * p[i - 1]
            P2 = ks * ms * (p[i - 1] + P1 * dt / 2)
            P3 = ks * ms * (p[i - 1] + P2 * dt / 2)
            P4 = ks * ms * (p[i - 1] + P3 * dt / 2)

    pp[i] = p[i] + (dt * (P1 + 2 * (P2 + P3) + P4) / 6)
    mms[i] = MMs * i

    return pp, mms


def block1(n):
    Xi = []
    Xs = []
    T = []
    P = []
    MMS = []

    # Начальные условия
    p0 = []
    for ii in range(n1):
        p0.append(0)
    P.append(p0)
    Xi.append(0)
    Xs.append(0)
    T.append(T0[n])

    j = 0

    while True:
        xi, xs, t, ms, J, ks, ki = kimpt1(Xi[j], Xs[j], T[j], n)
        p, mms = kfc1(ms, J, ks, ki, P[j])
        P.append(p)

        print(j, xs, xi, t, sep=" | ")

        Xs.append(xs)
        Xi.append(xi)
        T.append(t)
        MMS.append(mms)

        if Xs[j] >= 1:
            break
        j = j + 1

    return Xi, Xs, T, P, np.array(MMS)

block1(0)
