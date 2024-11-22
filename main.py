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

# Экспериментальные значения температуры
Texp1 = np.array(
    [
        [
            310.4,310.9,311.5,312.4,313.6,315.2,316.7,318.5,320.7,322.3,324.4,326.6,328.2,329.8,331.0,332.0,333.3,333.3,333.7,333.7,333.7,333.5,333.4,333.2,332.8,332.6,332.3,332.0,331.7,331.4,
        ],
        [309.6,310.0,310.3,310.7,311.4,312.0,313.0,314.0,315.0,316.2,317.5,319.0,320.4,322.3,323.7,325.4,326.8,327.9,329.0,330.0,330.7,331.2,331.6,331.8,332.0,332.0,331.9,331.9,331.6,331.4,
        ],
        [312.7,312.9,313.4,314.0,315.0,316.0,317.3,319.2,320.8,322.4,324.8,326.7,328.6,330.1,331.4,332.7,333.6,334.3,334.6,334.3,334.2,333.9,333.5,333.3,333.0,332.7,332.4,332.2,331.8,331.5,
        ],
        [309.3,309.4,309.6,310.0,310.3,310.9,311.6,312.2,313.3,314.2,315.7,317.0,318.3,319.8,321.6,323.1,324.5,326.3,327.4,328.2,329.2,329.9,330.4,330.8,331.1,331.1,331.1,331.0,330.9,330.7,
        ],
        [313.0,313.2,313.6,314.2,315.6,316.7,318.2,320.3,322.9,325.4,328.2,330.0,331.0,331.7,333.9,334.5,334.9,335.1,335.2,335.1,334.9,334.7,334.5,334.3,334.1,333.9,333.7,333.4,333.1,332.8,
        ],
        [311.2,311.4,311.8,312.6,313.4,314.4,315.5,316.6,317.7,319.4,321.0,322.4,324.2,325.7,327.4,328.8,330.1,331.0,332.0,332.5,333.1,333.5,333.7,333.8,333.9,333.8,333.7,333.5,333.5,333.1,
        ],
    ]
)


Texp2 = np.array(
    [
        [307.4,307.5,307.6,307.7,308.2,308.8,309.6,310.5,311.5,312.5,313.7,315.0,316.5,318.2,320.2,322.6,325.6,329.2,334.6,341.6,353.0,363.2,368.8,370.8,371.4,370.8,
        ],
        [310.4,310.6,311.3,312.2,313.3,314.4,315.9,317.4,319.4,321.4,323.9,327.3,332.1,338.6,349.6,361.4,368.9,372.7,373.9,373.9,373.5,372.6,371.8,
        ],
        [309.5,309.6,309.6,309.6,311.6,312.6,313.6,314.9,316.3,317.8,319.5,321.4,323.7,326.8,330.4,336.1,344.1,354.9,365.3,370.0,372.6,373.6,373.9,
        ],
        [313.3,314.3,315.4,316.7,318.5,320.5,323.2,325.8,329.5,335.9,343.8,353.9,361.3,367.4,371.0,373.1,373.5,373.8,373.6,373.1,372.6,
        ],
        [314.5,315.7,317.2,318.7,320.6,323.2,326.1,329.7,334.6,342.1,353.4,364.3,370.2,372.6,373.2,372.9,372.0,371.1,370.2,369.4,368.6,
        ],
        [309.4,309.6,309.8,310.2,310.8,312.1,313.1,314.3,315.7,317.1,318.9,321.2,323.7,326.9,331.0,337.9,346.8,356.4,363.0,368.1,371.2,372.7,373.3,
        ],
    ]
)


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

        Xs.append(xs)
        Xi.append(xi)
        T.append(t)
        MMS.append(mms)

        if Xs[j] >= 1:
            break
        j = j + 1

    return Xi, Xs, T, P, np.array(MMS)


def kimpt2(Xd, T, n):
    cps = Ms[n] / (Md2[n] + Ms[n] + Mr[n])
    cp = cps + ((Md2[n] * Xd) / (Mr[n]))
    Kf = A - (B * np.exp(C * cp))
    Kd = k * ((1 - b * cp) * (np.exp(-Ed / (R * T))))

    Xd1 = (Kd * np.power(Jac[n], 0.25 + wd * cp)) * (1 - Xd)
    T1 = (
        (Kted * Xd1 * md0[n] * Vt[n])
        - (
            ((Kf * Fst * L[n]) * (T - Th[n]) * Gh[n] * Ch * Ph)
            / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
        )
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    Xd2 = (Kd * np.power(Jac[n], 0.25 + wd * cp)) * (1 - (Xd + Xd1 * dt / 2))
    T2 = (
        (Kted * (Xd2) * md0[n] * Vt[n])
        - (
            (Kf * Fst * L[n] * ((T + T1 * dt / 2) - Th[n]) * Gh[n] * Ch * Ph)
            / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
        )
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    Xd3 = (Kd * np.power(Jac[n], 0.25 + wd * cp)) * (1 - (Xd + Xd2 * dt / 2))
    T3 = (
        (Kted * (Xd3) * md0[n] * Vt[n])
        - (
            (Kf * Fst * L[n] * ((T + T2 * dt / 2) - Th[n]) * Gh[n] * Ch * Ph)
            / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
        )
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    Xd4 = (Kd * np.power(Jac[n], 0.25 + wd * cp)) * (1 - (Xd + Xd3 * dt / 2))
    T4 = (
        (Kted * (Xd4) * md0[n] * Vt[n])
        - (
            (Kf * Fst * L[n] * ((T + T3 * dt / 2) - Th[n]) * Gh[n] * Ch * Ph)
            / (Kf * L[n] * Fst + Gh[n] * Ch * Ph)
        )
    ) / ((Map * Cap) + (Vt[n] * dpm * Cpm))

    xd = Xd + (dt * (Xd1 + 2 * (Xd2 + Xd3) + Xd4) / 6)
    t = T + (dt * (T1 + 2 * (T2 + T3) + T4) / 6)
    md = md0[n] * (1 - xd)
    kd = Kd / (np.power(Jac[n], (0.75 - (wd * cp))))

    return md, kd, xd, t


def kfc2(kd, md, p):
    pp = np.zeros(n1)
    mmn = np.zeros(n1)

    for i in range(n1):
        if i == 0:
            P1 = -kd * md * p[i]
            P2 = -kd * md * (p[i])
            P3 = -kd * md * (p[i])
            P4 = -kd * md * (p[i])

        elif i > 0 and i < n1 - 2:
            P1 = -kd * md * (p[i]) + kd * md * p[i - 1]
            P2 = -kd * md * (p[i]) + kd * md * p[i - 1]
            P3 = -kd * md * (p[i]) + kd * md * p[i - 1]
            P4 = -kd * md * (p[i]) + kd * md * p[i - 1]

        else:
            P1 = kd * md * p[n1 - 1]
            P2 = kd * md * p[n1 - 1]
            P3 = kd * md * p[n1 - 1]
            P4 = kd * md * p[n1 - 1]

        pp[i] = p[i] + (dt * (P1 + 2 * (P2 + P3) + P4) / 6)
        if pp[i] < 0:
            pp[i] = 0
        mmn[i] = ns * MMs + (i - ns) * MMd

    return pp, mmn


def block2(n):
    Xd = []
    T2 = []
    P2 = []
    MM = []

    p2 = P[N - 1]
    P2.append(p2)

    # Начальные условия
    Xd.append(0)
    T2.append(Td0[n])

    j = 0

    while True:
        md, kd, xd, t = kimpt2(Xd[j], T2[j], n)
        pp, mm = kfc2(kd, md, P2[j])
        Xd.append(xd)
        T2.append(t)
        P2.append(pp)
        MM.append(mm)

        if Xd[j] >= 1:
            break
        else:
            j = j + 1

    return np.array(Xd), np.array(T2), np.array(P2), np.array(MM)


Xd, T2, P2, MM = block2(1)


cwd = os.getcwd()


N = len(P2)


y = np.arange(0, 30, 1)
x = np.arange(0, N, 1)
X, Y = np.meshgrid(x, y)


# Pn = pd.DataFrame(P2[N-1])
# Mn = pd.DataFrame(MM[N-2])
# Pn.to_csv(os.path.join(cwd,"P2.csv"))
# Mn.to_csv(os.path.join(cwd,"M2.csv"))


plt.plot(Xd, label="Конверсия мономера")
plt.xlabel("Время, мин")
plt.ylabel("Конверсия, ед")
plt.legend()
plt.show()
# plt.plot(P2[N-1][:50])

nn = len(T2)


plt.plot(T2, label="Расчетные значение температура")
plt.plot(Texp2[1][:nn], "b*", label="Экспериментальные данные")
plt.xlabel("время, мин")
plt.ylabel("Температура, К")
plt.legend()
plt.show()


# plt.plot(P2[N-1,:100])
# plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.plot_surface(X, Y, P2[:, 5:35].T)
plt.xlabel("номер фракции")
plt.ylabel("Время, мин")
ax.set_zlabel("концентрации,      моль/л")


plt.show()


import pandas as pd
import matplotlib.pyplot as plt
from block1 import *
import os


Xi, Xs, T, P, MMS = block1(0)
P = np.array(P)
# PP = pd.DataFrame(P)


N = len(P)
T = np.array(T)
# Plot the 3d graph for P


y = np.arange(0, 30, 1)
x = np.arange(0, N, 1)
X, Y = np.meshgrid(x, y)


cwd = os.getcwd()
path = cwd


# Plot the concentration
# N = len(P)


# for i in range(len(P)-2):
# plt.plot(MMS[i][:20], P[i][:20])
# plt.show()
# for i in range(len(P)-10,len(P)-2):
# plt.plot( MMS[i][:20],P[i][:20])
# plt.show()
# plt.plot(MMS[N-2][:30], P[N-1][:30])
# plt.show()
# Pn = pd.DataFrame(P)
# Pn.to_csv(os.path.join(cwd,"csv_data1.csv"))
# for i in range(N):
# plt.plot(P[i])


# P[i][j]
# Конверсия
# plt.plot(Xi, label="Xi -конверсия катализатора")
# plt.plot(Xs, label="Xs -конверсия мономера")
# plt.xlabel("Время")
# plt.ylabel("Конверсия")
# plt.legend()
# plt.show()


nn = len(Texp1[0])
plt.plot(T[:nn], label="Расчетные значение температура")
plt.plot(Texp1[0], "b*", label="Экспериментальные данные")
plt.xlabel("время, мин")
plt.ylabel("Температура, K")
plt.legend()
plt.show()

nn = len(Texp1[0])
plt.plot(Xi, label="Xi")
plt.plot(Xs, label="Xs")
plt.xlabel("Время, мин")
plt.ylabel("конверсия, ед")
plt.legend()
plt.show()

# plt.plot(P[N-1][:50])
# plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")


ax.plot_surface(X, Y, P[:, 5:35].T)
plt.xlabel("номер фракции")
plt.ylabel("Время, мин")
ax.set_zlabel("концентрации, моль/л", labelpad=15.0)


plt.show()


from graph2 import P2, MM, N


np.random.seed(42)
P = P2[N - 1]
M = MM[N - 2]


# Расчет матрица A
def A(P, MM, N, dP):

    a = []
    for n in range(N):

        mn = int(P[n] / dP)
        for i in range(mn):
            a.append(MM[n])

    return np.array(a)


# Расчет матрица B
np.random.seed(42)


def B(a):

    # S = 0
    N = len(a)
    b = []
    for j in range(N):
        S = 0
        for k in range(4):
            i = np.random.randint(N)
            S += S + a[i]
            np.delete(a, i)
        b.append(S)

        if len(a) < 4:
            break
    return np.array(b)


# Расчет матрица C
np.random.seed(42)


def C(b, MM, dP, dMM, N):
    C = np.zeros((N, 2))
    for i in range(N):
        S = 0
        NN = 1

    for j in range(len(b)):
        if (b[j] > MM[i]) and (b[j] <= MM[i] + dMM):
            S = S + b[j]
            NN += 1
            # print(S)

    C[i][0] = S / NN
    C[i][1] = NN * dP

    return np.array(C)


# Расчет соедные значение молекулрные массы и концентрации
def final_C(C, K):
    c = np.zeros((len(C), 2))
    for k in range(K):
        for i in range(len(C) - 1):
            c[i][0] = (C[i][0] + C[i + 1][0]) / 2
            c[i][1] = (C[i][1] + C[i + 1][1]) / 2

    return np.array(c)


# Расчет показатели
def final_values(P, MM):

    sum_mn = sum([(MM[n] * P[n]) for n in range(len(P))])
    sum_p = sum([P[n] for n in range(len(P))])
    sum_mn2 = sum([((MM[n] ** 2) * P[n]) for n in range(len(P))])

    mn = sum_mn / sum_p
    mw = sum_mn2 / sum_mn

    return mn, mw


# Основная функция


def pokasateli(M, P):
    N = len(P)
    n = int(N / 2)

    # sum of P
    sum_p = sum([P[i] for i in range(N)])
    sum_mp = sum([M[i] * P[i] for i in range(N)])
    sum_m2p = sum([(M[i] ** 2) * P[i] for i in range(N)])

    mn = sum_mp / sum_p
    mw = sum_m2p / sum_mp

    sum_mwp = sum([((M[i] * P[i]) - (mn * P[n])) ** 2 for i in range(N)])
    sigma = np.sqrt((1 / (N - 1)) * sum_mwp)
    kd = mw / mn

    return mn, mw, sigma, kd


def MMP(P, MM, dP, dMM, K=4):
    N = len(P)

    a = A(P, MM, N, dP)
    b = B(a)
    c = C(b, MM, dP, dMM, N)
    CC = final_C(c, K)
    mn, mw = pokasateli(CC[:, 0], CC[:, 1])

    return mn, mw, CC


# Вычичления
_, _, CMMP = MMP(P, M, dP, dM)
plt.subplot(1, 2, 1)
plt.plot(CMMP[550:750, 0], CMMP[550:750, 1])
plt.title("результат имитационного модели ММР")


# plt.show()


plt.subplot(1, 2, 2)
plt.plot(M[5:35], P[5:35], "r")  #
plt.title("Рассчитанный модель")
plt.show()
