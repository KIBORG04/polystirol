from scipy.stats import poisson, norm
import numpy as np

# Параметр λ
lam = 250

# Диапазон значений
x = np.arange(0, 1000)  # Замените 20 на нужное максимальное значение

# Расчёт вероятностей
probabilities = poisson.pmf(x, mu=lam)

mu = 0
sigma = 100

# Расчёт вероятностей
#probabilities = norma.pdf(x, loc=mu, scale=sigma)

print(probabilities)