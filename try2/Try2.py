from constants import *
from my_time import *
import block1, block2, my_time
import math
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import scipy.integrate as sci

# Номер опыта
n = 0

def main():
    values = block1.block1(n)
    printBlock1Stats(values)

    P = [ar[-1] for ar in values.y[3:]]
    values2 = block2.block2(n, P)
    printBlock2Stats(values2)

    printTemperatureChart(values, values2)

def printTemperatureChart(sol1, sol2):
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=["Сравнение для первой модели", "Сравнение для второй модели"]
    )
    fig.add_trace(go.Scatter(x=np.arange(timeStart, timeEnd, 1), y=Texp1[n], mode='markers', name='Экспериментальная температура',
                            marker=dict(color='blue', size=8)),
                            row=1, col=1)
    fig.add_trace(go.Scatter(x=sol1.t, y=sol1.y[2], mode='lines', name='Температура модели',
                            line=dict(color='red', width=2)),
                            row=1, col=1)
    fig.add_trace(go.Scatter(x=np.arange(timeStart, timeEnd, 2), y=Texp2[n], mode='markers', name='Экспериментальная температура',
                            marker=dict(color='blue', size=8)),
                            row=2, col=1)
    fig.add_trace(go.Scatter(x=sol2.t, y=sol2.y[1], mode='lines', name='Температура модели',
                            line=dict(color='red', width=2)),
                            row=2, col=1)
    # Общие настройки оформления
    fig.update_layout(
        height=1600,  # Высота графика
        title="Сравнение экспериментальной и смоделированной температуры",
        xaxis_title="Time (s)",
        yaxis_title="Temperature (°C)",
        legend=dict(title="Data"),
        template="plotly_white"
    )
    fig.show()

def printBlock1Stats(sol):
    print(f'Конверсии: Xi:{sol.y[0][-1]} Xs:{sol.y[1][-1]}')
    print(f'Отклонение температуры на стадии 1: {sum([(sol.y[2][int(i*1/dt)]-Texp1[n][i]) for i in range(len(Texp1[n]))])}')
    print([ar[-1] for ar in sol.y[3:]])

def printBlock2Stats(sol):
    print(f'Конверсии: Xd:{sol.y[0][-1]}')
    print(f'Отклонение температуры на стадии 2: {sum([(sol.y[1][int(i*2/dt)]-Texp2[n][i]) for i in range(len(Texp2[n]))])}')
    print([ar[-1] for ar in sol.y[2:]])

if __name__ == "__main__":
    main()