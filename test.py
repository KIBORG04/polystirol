import numpy as np

# Исходные параметры
F = 5
E1 = 0.5
D = 1
E2 = 0.45
R = 4
epsilon = 0.005
k = 2
T0 = 105 + 273.15
NF = int(16 * E2)
N = int(80 * E2)
TN = 75 + 273.15
P0 = 1.25

XF = np.array([20, 80])

C = np.array([
    [123.912, 70.435],
    [-8753.9, -7362.7],
    [0, 0],
    [0.020198, 0.006952],
    [0, 0],
    [-18.1, -9]
])

def getL(R):
    L = np.zeros((N+2))
    for i in range(N+2):
        if i == 0:
            L[i] = F-D
        elif 1 <= i and i <= NF:
            L[i] = R+F
        elif NF < i and i <= N+1:
            L[i] = R
    return L

def getV(R):
    V = np.zeros((N+2))
    for i in range(N+2):
        V[i] = R+D
    return V

# Функция для вычисления температуры
def calculate_T(T0):
    return np.array([T0 - ((T0 - TN) * i / N) for i in range(N+1)])

# Функция для вычисления Y
def calculate_Y(T):
    a12 = np.exp(643.7 / T - 143350 / T**2)
    b12 = np.exp(112.4714 / T - 264382.77 / T**2)
    a21 = 1 / a12
    b21 = b12 / a12
    c = 0.2732

    y0 = 10**(c * b12 * XF[1]**2 * (a12 * XF[1] + (1 + XF[0]) + (b12-1)*XF[1]**2)/(XF[0]**2+a12*XF[1]**2+b12*XF[0]*XF[1])**2)
    y1 = 10**(c * b21 * XF[0]**2 * ((a21 * XF[0]) + (1 + XF[1]) + (b21-1)* XF[1]**2)/(XF[1]**2+a21*XF[0]**2+b21*XF[0]*XF[1])**2)

    return np.array([y0, y1])

# Функция для вычисления P
def calculate_P(T):
    P = np.zeros((N+1, 2))
    for i in range(N+1):
        for j in range(2):
            P[i, j] = np.exp(
                C[0, j] + C[1, j] / (C[2,j]+T) + C[3, j] * T + C[4, j] * T**2 + C[5, j] * np.log(T)
            )
    return P

# Функция для вычисления α
def calculate_alpha(T):
    alpha = np.zeros((N+1, 2))
    for i in range(N+1):
        Y = calculate_Y(T[i])
        for j in range(2):
            P = calculate_P(T[i])
            alpha[i, j] = Y[j] * P[i, j]
    return alpha

# Функция для вычисления K
def calculate_K(X, T0):
    alpha = calculate_alpha(calculate_T(T0))
    K = np.zeros((N+1, 2))
    for i in range(N+1):
        for j in range(2):
            K[i, j] = alpha[i, j] / np.sum(alpha[i, :] * X[i, :])
    return K

# Нормализация матрицы X
def normalize_X(X):
    Xnorm = np.zeros_like(X)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            Xnorm[i, j] = X[i, j] / np.sum(X[i, :])
    return Xnorm

def calculate_matrix_X(K, R, D, XF):
    L = getL(R)
    V = getV(R)
    X = np.zeros((N + 2, 2))
    G = np.zeros((N + 2, N + 2))

    # Построение матрицы G и массива d
    a = np.zeros(N+2)
    b = np.zeros(N+2)
    c = np.zeros(N+2)
    d = np.zeros(N+2)
    for j in range(2):
        for i in range(N + 2):
            if i > 0:
                a[i] = V[i-1]*K[i-1, j]
            if i < N + 1:
                b[i] = -(L[i]+V[i]*K[i,j])
                c[i]=L[i+1]
            if i == NF:
                d[i] = -F * XF[j]
            if i == N+1:
                b[i] = -(R+D)
                c[i] = 0
        for i in range(N + 2):
            if i > 0:
                G[i, i-1] = a[i]
            G[i,i] = b[i]
            if i < N+1:
                G[i,i+1] = c[i]
        X[:,j] = np.linalg.inv(G) @ d

    return X


# Основной алгоритм Program
async def program(state, R, D, T0):
    n = 0
    g = 0
    z = 0
    X1 = np.tile(XF, (N + 1, 1))

    while g == 0:
        n += 1
        K1 = calculate_K(X1,T0)
        Xres = calculate_matrix_X(K1, R, D, XF)
        g = check_convergence(Xres, z)
        X1 = normalize_X(Xres)
        z += 1

    return X1

# Функция для проверки сходимости
def check_convergence(X, z):
    sum = 0
    for row in X:
        if np.abs(1 - np.sum(row)) < epsilon or z == 13:
            sum += 1
        else:
            return 0
    return len(X) == sum

import asyncio

async def program2222(RR, T00, max_iter=50, A=0.5, B=0.5, xi=0.1):
    S = np.zeros(max_iter)
    S2 = np.zeros(max_iter)
    S3 = np.zeros(max_iter)
    dS = np.zeros(max_iter)
    dS2 = np.zeros(max_iter)
    Q = np.zeros(3)
    Kr = np.zeros(max_iter)

    for i in range(max_iter):
        # Вычисление X на основе текущих значений RR и T00
        X_task = program("x", RR, D, T00)
        X2_task = program("x", RR + xi * RR, D, T00)
        X3_task = program("x", RR, D, T00 + xi * T00)

        X, X2, X3 = await asyncio.gather(X_task, X2_task, X3_task)

        S[i] = A * X[0][0] + B * X[N + 1][1]
        S2[i] = A * X2[0][0] + B * X2[N + 1][1]
        S3[i] = A * X3[0][0] + B * X3[N + 1][1]

        dS[i] = (S2[i] - S[i]) / (xi * S[i])
        dS2[i] = (S3[i] - S[i]) / (xi * S[i])

        # Вычисление d
        d = np.sqrt(dS[i]**2 + dS2[i]**2)

        # Условие выхода
        if d < 0.01:
            Q[0] = RR
            Q[1] = T00
            Q[2] = i
            break

        RR = RR - xi * (dS[i] / d)
        T00 = T00 - xi * (dS2[i] / d)

        Q[0] = RR
        Q[1] = T00
        Q[2] = i

        Kr[i] = S[i]

    # Новый вызов программы
    newX = await asyncio.gather(program("x", RR, D, T00))

    return Q, newX[0]

async def getbased():
    newX = await asyncio.gather(program("x", R, D, T0))
    return newX[0]

import plotly.graph_objects as go
# Программа
result = asyncio.run(program2222(R, T0))
newX = asyncio.run(getbased())
print(f'Коэффициенты до: R: {R} T0: {T0}')
print(f'Коэффициенты после: R: {result[0][0]} T0: {result[0][1]}')
print(f'Voda xN: {newX[-1][1]} Spirt x0: {newX[0][0]}')
print('Оптимальное')
print(f'Voda xN: {result[1][-1][1]} Spirt x0: {result[1][0][0]}')
result = result[1]
#print("Результат программы:", result)
#print(newX)
fig = go.Figure()
fig.add_trace(go.Scatter(x=np.arange(0, len(newX)), y=newX[:, 0], name='Вода'))
fig.add_trace(go.Scatter(x=np.arange(0, len(newX)), y=newX[:, 1], name='Спирт', mode="lines+markers"))
fig.show()
fig = go.Figure()
fig.add_trace(go.Scatter(x=np.arange(0, len(result)), y=result[:, 0], name='Вода'))
fig.add_trace(go.Scatter(x=np.arange(0, len(result)), y=result[:, 1], name='Спирт', mode="lines+markers"))
fig.show()