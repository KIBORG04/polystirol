import numpy as np

# Текущее время
timeStart = 0
# Время моделирвания
timeEnd = 75
# Шаг дискретизации времени
TIMESTEP = 0.01
# Короткий псевдоним
dt = TIMESTEP

# Все дискретные шаги времени
time = np.arange(timeStart, timeEnd, TIMESTEP)
# Количество шагов
steps = len(time)