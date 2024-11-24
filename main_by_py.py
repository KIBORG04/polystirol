import math
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

# Количество фракций
FACTIONS = 6

# Начальная температура для каждой фракции
T0 = [310.4, 309.6, 312.7, 309.3, 313.0, 311.2]
# Начальная концентрация стирола, моль/л.
ms0 = [0.745, 0.7304, 0.7129, 0.7106, 0.6683, 0.7168]
# Концентрация катализатора , моль/л.
Jk = [0.00478, 0.0045, 0.00416, 0.00361, 0.0041, 0.00402]
# Масса растворителя мономеров
Ms = [1140, 1080, 1085, 1081, 998, 1085]
# Масса растворителя мономеров
Mr = [11400, 11140, 11370, 11500, 11270, 11290]
# Температура хладагента Тh, °С
Th = [290, 289.5, 288, 290, 296, 288.5]
# Расход хладагента Gh, л/мин
Gh = [1500, 2000, 2000, 1666, 1500, 1500]
# Молекулярная масса стирола, кг/моль
mms = 0.092
MMs = [mms * i for i in range(FACTIONS)]
# Молекулярная масса дивинила, кг/моль
MMd = 0.054

# Динамическая вязкость
ng = [0.125, 0.129, 0.13, 0.147, 0.129, 0.137]
Mr = [11400, 11140, 11370, 11500, 11270, 11290]
Md1 =[1765, 1687, 1701, 1664, 1558, 1702]
Md2 =[837, 715, 715, 722, 663, 710]
pr = [0.748, 0.748, 0.739, 0.75, 0.744, 0.739]
ps = [0.91, 0.916, 0.92, 0.919, 0.921, 0.92]
pd = [0.638, 0.637, 0.637, 0.635, 0.637, 0.637]
Td0 =[307.4, 310.4, 309.5, 313.3, 314.5, 309.4]

# Что-то
Vk = 139.2
# объём реакционной смеси, l
Vt = [(Vk + ((Ms[i] / ps[i]) + (Mr[i] / pr[i]))) for i in range(len(Ms))]
# расчет степени заполнения
L = [((Vt[i] - 3000) / 28500) for i in range(len(Vt))]

Vd = np.array([(Md1[i] / pd[i]) for i in range(len(pd))])
# Концентрация активных центровполимеризации на первой и второй стадии синтеза
Jac = np.array([(Jk[i] * ((Vt[i]) / (Vt[i] + Vd[i]))) for i in range(len(Vt))])
md0 = np.array([((pd[i] * Vd[i]) / (MMd * (Vt[i] + Vd[i]))) for i in range(len(pd))])

# Постоянные констант скорости инициирования и роста цепи [л/(мин*моль)]
Ki0, Ks0, Kd0 = 0.835 * 10**10, 5.76 * 10**11, 2 * 10**12
# Энергия активации инициирования и роста полистирольных цепей [Дж/моль]
Ei, Es, Ed = 59962, 71184, 80850
# Универсальная газовая постоянная
#R = 8.31446
R = 8.32
#
b = 3.77

# Параметр зависимости, учитывающей влияние концентрации полимера на порядок реакции по катализатору
ws = -0.4142
#
wd = -0.184
# Масса аппарата (реактора), кг
Map = 16840
# Теплоёмкость материала аппарата, Дж/(кг*К)
Cap = 0.46
# Плотность реакционной массы, кг/л
dpm = 0.763
# Теплоёмкость реакционной массы, кДж/(кг*К)
Cpm = 1.8436

# Теплоёмкость (кДж/(кг*К)) и плотность хладагента (кг/л)
Ch, ph = 4.19, 1

# Коэффициенты идентификации [Дж/(м2*К*мин)]
A, B, C = 14.17, 1.007, 10.38

# Площадь поверхности теплообмена в реакторе, m^2
Fst = 42.0
# Коэффициент тепловедение стирола и дивинила
Kted, Ktes = 80.9, 74.87

# Экспериментальные значения температуры
Texp1 = [
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

Texp2 = [
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

# Текущее время
timeStart = 0.0
# Время моделирвания
timeEnd = 35
# Шаг дискретизации времени
TIMESTEP = 0.01
# Короткий псевдоним
dt = TIMESTEP

# Все дискретные шаги времени
time = np.arange(timeStart, timeEnd, TIMESTEP)
# Количество шагов
steps = len(time)


def kimpt(prevXi, prevXs, prevT, n):
    Cp = getCp(prevXs, n)

    # Коэффициент теплопередачи через стенку аппарата
    Kf = A - B * math.pow(math.e, C * Cp)

    Ki = getKi(prevT)
    Ks = getKs(prevT)

    def getNewXi(xi, xs):
        return Ki * ms0[n] * (1 - xi) * (1 - xs)

    def getNewXs(xi, xs):
        return (Ki * Jk[n] * (1 - xi) + Ks * ((Jk[n] * xi)**(0.5 + ws*Cp)))*(1-xs)

    # Температура
    def getNewT(xs, t):
        numeratorLeft = (Ktes * xs * ms0[n] * Vt[n])
        numeratorRigth = (Kf * Fst * L[n] * (t - Th[n]) * Gh[n] * Ch * ph)/(Kf * L[n] * Fst + Gh[n] * Ch * ph)
        denominator = (Map * Cap) + (Vt[n] * dpm * Cpm)
        return (numeratorLeft - numeratorRigth) / denominator

    # Попробовать заменить dXs1 и другие d-шки в расчетах температуры и передавать туда уже посчитанное newXs

    dXi1 = dt * getNewXi(prevXi, prevXs)
    dXs1 = dt * getNewXs(prevXi, prevXs)

    dXi2 = dt * getNewXi(prevXi + dXi1/2, prevXs + dXs1/2)
    dXs2 = dt * getNewXs(prevXi + dXi1/2, prevXs + dXs1/2)

    dXi3 = dt * getNewXi(prevXi + dXi2/2, prevXs + dXs2/2)
    dXs3 = dt * getNewXs(prevXi + dXi2/2, prevXs + dXs2/2)

    dXi4 = dt * getNewXi(prevXi + dXi3, prevXs + dXs3)
    dXs4 = dt * getNewXs(prevXi + dXi3, prevXs + dXs3)

    newXi = prevXi + (dXi1 + 2*dXi2 + 2*dXi3 + dXi4)/6
    newXs = prevXs + (dXs1 + 2*dXs2 + 2*dXs3 + dXs4)/6

    dT1  = dt * getNewT(dXs1/dt, prevT)
    dT2  = dt * getNewT(dXs2/dt, prevT + dT1/2)
    dT3  = dt * getNewT(dXs3/dt, prevT + dT2/2)
    dT4  = dt * getNewT(dXs4/dt, prevT + dT3)

    #print(f"Xi1: {dXi1} Xi2: {dXi2} Xi3: {dXi3} Xi4: {dXi4}")
    #print(f"Xs1: {dXs1} Xi2: {dXs2} Xs3: {dXs3} Xs4: {dXs4}")
    #print(f"T1: {dT1} T2: {dT2} T3: {dT3} T4: {dT4}")

    newT  = prevT +  (dT1 +  2*dT2  + 2*dT3  + dT4) /6

    return newXi, newXs, newT

def kfcModel(prevP, xs, xi, t, n):
    def getP0(xi, ms, ki, ks, p):
        J = Jk[0] * (1 - xi)
        return -ks * ms * p + ki * ms * J

    def getPn(ms, ks, pn, prevPn):
        return -ks * ms * pn + ks * ms * prevPn

    def getPN(ms, ks, p):
        return getPn(ms, ks, 0, p)

    n = n + 1
    newP = [0.0 for _ in range(n)]
    for i in range(n):
        ms = ms0[i] * (1 - xs)
        Cp = getCp(xs, i)
        ks = getKs(t)/(Jk[i]**(0.5-ws*Cp))
        ki = getKi(t)
        if i == 0:
            dP1 = dt * getP0(xi, ms, ki, ks, prevP[0])
            dP2 = dt * getP0(xi, ms, ki, ks, prevP[0] + dP1/2)
            dP3 = dt * getP0(xi, ms, ki, ks, prevP[0] + dP2/2)
            dP4 = dt * getP0(xi, ms, ki, ks, prevP[0] + dP3)
        elif i == n - 1:
            dP1 = dt * getPN(ms, ks, prevP[i-1])
            dP2 = dt * getPN(ms, ks, prevP[i-1] + dP1/2)
            dP3 = dt * getPN(ms, ks, prevP[i-1] + dP2/2)
            dP4 = dt * getPN(ms, ks, prevP[i-1] + dP3)
        else:
            dP1 = dt * getPn(ms, ks, prevP[i],         prevP[i-1])
            dP2 = dt * getPn(ms, ks, prevP[i] + dP1/2, prevP[i-1] + dP1/2)
            dP3 = dt * getPn(ms, ks, prevP[i] + dP2/2, prevP[i-1] + dP2/2)
            dP4 = dt * getPn(ms, ks, prevP[i] + dP3  , prevP[i-1] + dP3)
        newP[i] = prevP[i] + (dP1 + 2*dP2 + 2*dP3 + dP4)/6
    return newP

def getCp(xs, n):
    return (Ms[n] * xs) / Mr[n]

# Константа скорости инициирования
def getKi(t):
    return Ki0 * math.pow(math.e, -Ei/(R*t))

# Константа скорости роста
def getKs(t):
    return Ks0 * math.pow(math.e, -Es/(R*t))

def block1(n):
    # Фракций не должно быть 0, но индексы листов начинаются с 0
    n = n - 1

    # Конверсия катализатора во времени
    Xi = [0]
    # Конверсия мономера (стирола) во времени
    Xs = [0]
    # Температура во времени
    T = [T0[n]]
    # Концентрация молекул каждой фракции
    P = []

    pp = [0 for _ in range(n+1)]
    P.append(pp)

    for i in range(steps):
        xi, xs, t = kimpt(Xi[i], Xs[i], T[i], n)
        p = kfcModel(P[i], xs, xi, t, n)

        Xi.append(xi)
        Xs.append(xs)
        T.append(t)
        P.append(p)

        #print(xi, xs, t, sep="\t| ")
        #print(p)
        if xs >= 1:
            break

    time_exp = np.arange(timeStart, timeEnd, 1)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=time_exp, y=Texp1[n], mode='markers', name='Экспериментальная температура',
                            marker=dict(color='blue', size=8)))
    fig.add_trace(go.Scatter(x=time, y=T, mode='lines', name='Температура модели',
                            line=dict(color='red', width=2)))
    fig.update_layout(
        title="Сравнение экспериментальной и смоделированной температуры",
        xaxis_title="Time (s)",
        yaxis_title="Temperature (°C)",
        legend=dict(title="Data"),
        template="plotly_white"
    )
    fig.show()
    return P

def kimpt2(prevXd, prevT, n):
    cps = Ms[n] / (Md2[n] + Ms[n] + Mr[n])
    Cp = cps + (Md2[n] * prevXd) / Mr[n]

    Kf = A - B * math.pow(math.e, C * Cp)
    Kd = Kd0 * (1 - b * Cp) * math.pow(math.e, -Ed / (R * prevT))

    def Xd(xd):
        return (Kd * math.pow(Jac[n], 0.25 + wd*Cp)) * (1 - xd)

    def T(t, xd):
        numeratorLeft = (Kted * xd * md0[n] * Vt[n])
        numeratorRigth = (Kf * Fst * L[n] * (t - Th[n]) * Gh[n] * Ch * ph)/(Kf * L[n] * Fst + Gh[n] * Ch * ph)
        denominator = (Map * Cap) + (Vt[n] * dpm * Cpm)
        return (numeratorLeft - numeratorRigth) / denominator

    dXd1 = dt * Xd(prevXd)
    dXd2 = dt * Xd(prevXd + dXd1/2)
    dXd3 = dt * Xd(prevXd + dXd2/2)
    dXd4 = dt * Xd(prevXd + dXd3)
    newXd = prevXd + (dXd1 + 2*dXd2 + 2*dXd3 + dXd4) / 6

    dT1 = dt * T(prevT        , dXd1/dt)
    dT2 = dt * T(prevT + dT1/2, dXd2/dt)
    dT3 = dt * T(prevT + dT2/2, dXd3/dt)
    dT4 = dt * T(prevT + dT3  , dXd4/dt)
    newT  = prevT + (dT1 + 2*dT2 + 2*dT3 + dT4) / 6

    return newXd, newT

def kfcModel2(prevP, xd, t, n):
    def getP0(kd, md, p):
        return -kd * md * p

    def getPn(kd, md, pn, prevPn):
        return -kd * md * pn + kd * md * prevPn

    def getPN(kd, md, p):
        return getPn(kd, md, 0, p)

    cps = Ms[n] / (Md2[n] + Ms[n] + Mr[n])
    Cp = cps + (Md2[n] * xd) / Mr[n]

    n = n + 1
    newP = [0.0 for _ in range(n)]
    for i in range(n):
        md = md0[i]*(1-xd)
        # TODO: я в формулу сам подставил (1 - b * Cp) !!! В модели этого нету. Результат: снижение концентрации молекул в много раз
        kd = (Kd0 * (1 - b * Cp) * math.pow(math.e, -Ed/(R*t))) / (Jac[i]**(0.75 - wd*Cp))
        if i == 0:
            dP1 = dt * getP0(kd, md, prevP[0])
            dP2 = dt * getP0(kd, md, prevP[0] + dP1/2)
            dP3 = dt * getP0(kd, md, prevP[0] + dP2/2)
            dP4 = dt * getP0(kd, md, prevP[0] + dP3)
        elif i == n - 1:
            dP1 = dt * getPN(kd, md, prevP[i-1])
            dP2 = dt * getPN(kd, md, prevP[i-1] + dP1/2)
            dP3 = dt * getPN(kd, md, prevP[i-1] + dP2/2)
            dP4 = dt * getPN(kd, md, prevP[i-1] + dP3)
        else:
            dP1 = dt * getPn(kd, md, prevP[i],         prevP[i-1])
            dP2 = dt * getPn(kd, md, prevP[i] + dP1/2, prevP[i-1] + dP1/2)
            dP3 = dt * getPn(kd, md, prevP[i] + dP2/2, prevP[i-1] + dP2/2)
            dP4 = dt * getPn(kd, md, prevP[i] + dP3,   prevP[i-1] + dP3)
        newP[i] = prevP[i] + (dP1 + 2*dP2 + 2*dP3 + dP4)/6
    return newP

def block2(P1, n):
    # Фракций не должно быть 0, но индексы листов начинаются с 0
    n = n - 1

    Xd = [0]
    T = [T0[n]]
    # Концентрация молекул каждой фракции
    P = [P1[len(P1)-1]]
    # Моли
    MMn = []

    for i in range(steps):
        xd, t = kimpt2(Xd[i], T[i], n)
        p = kfcModel2(P[i], xd, t, n)

        Xd.append(xd)
        T.append(t)
        P.append(p)

        #print(xd, t, p, sep=" | ")
        #print(p)
        if xd >= 1:
            break

    time_exp = np.arange(timeStart, timeEnd, 1)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=time_exp, y=Texp2[n], mode='markers', name='Экспериментальная температура',
                            marker=dict(color='blue', size=8)))
    fig.add_trace(go.Scatter(x=time, y=T, mode='lines', name='Температура модели',
                            line=dict(color='red', width=2)))
    fig.update_layout(
        title="Сравнение экспериментальной и смоделированной температуры",
        xaxis_title="Time (s)",
        yaxis_title="Temperature (°C)",
        legend=dict(title="Data"),
        template="plotly_white"
    )
    fig.show()

    fig = go.Figure()

    for i in range(n+1):
        fig.add_trace(go.Scatter(
            x=time, 
            y=[row[i] for row in P], 
            mode='lines', 
            name=f'Фракция {i + 1}',
            marker=dict(color='blue', size=8)
        ))

    fig.update_layout(
        title="Концентрация",
        xaxis_title="Time (s)",
        yaxis_title="Моль на л",
        legend=dict(title="Data"),
        template="plotly_white"
    )
    #fig.show()

    return P, MMn

# Среднечисленная молекулярная масса:
def getMn(P, MMn):
    numerator = 0
    denumerator = 0
    for i in range(FACTIONS):
        numerator += P[i] * MMn[i]
        denumerator += P[i]
    return numerator / denumerator

# Средневзвешенная молекулярная масса:
def getMv(P, MMn):
    numerator = 0
    denumerator = 0
    for i in range(FACTIONS):
        numerator += P[i] * MMn[i]**2
        denumerator += P[i] * MMn[i]
    return numerator / denumerator

# Среднеквадратическое отклонение ММР
def getStdVar(P, MMn):
    sum = 0
    Mn = getMn(P, MMn)
    for i in range(FACTIONS):
        sum = (MMs[i] * P[i] - Mn * P[FACTIONS / 2])**2 # почему FACTIONS / 2 ?
    return math.sqrt(sum/(FACTIONS - 1))

def getA(P, MM, dP):
    a = []
    for n in range(len(P)):
        mn = int(P[n] / dP)
        for _ in range(mn):
            a.append(MM[n])
    return np.array(a)

def getB(A):
    N = len(A)
    b = []
    for _ in range(N):
        S = 0
        for _ in range(4):
            i = np.random.randint(N)
            S += S + A[i]
            np.delete(A, i)
        b.append(S)

        if len(A) < 4:
            break
    return np.array(b)

def getC(B, MM, N, dP, dMM):
    C = np.zeros((N, 2))
    for i in range(N):
        sum = 0
        NN = 1
        for j in range(len(B)):
            if (B[j] > MM[i]) and (B[j] <= MM[i] + dMM):
                sum = sum + B[j]
                NN += 1

        C[i][0] = sum / NN
        C[i][1] = NN * dP

    # общее количество повторений процедуры усреднения
    K = 4
    c = np.zeros((len(C), 2))
    for k in range(K):
        for i in range(len(C) - 1):
            c[i][0] = (C[i][0] + C[i + 1][0]) / 2
            c[i][1] = (C[i][1] + C[i + 1][1]) / 2
    return C

def calcMatrixC(P, MMn, dP=3, dMM=3):
    A = getA(P, MMn, dP)
    B = getB(A)
    C = getC(B, MMn, len(P), dP, dMM)
    return C

def printCharacteristics(C):
    pp = C[:, 0]
    mmn = C[:, 1]
    Mn = getMn(pp, mmn)
    Mv = getMv(pp, mmn)
    sigma = getStdVar(pp, mmn)

    print(f"Mn: {Mn}; Mv: {Mv}; sigma: {sigma}")
    return C

n = 6
P = block1(n)
P2, MMn = block2(P, n)

#C = calcMatrixC(P2[len(P2)-1], MMn[len(MMn)-1])
#printCharacteristics(C)
#
#plt.subplot(1, 2, 1)
#plt.plot(C[550:750, 0], C[550:750, 1])
#plt.title("результат имитационного модели ММР")
#
#plt.subplot(1, 2, 2)
#plt.plot(MMn[5:35], P2[5:35], "r")  #
#plt.title("Рассчитанный модель")
#plt.show()