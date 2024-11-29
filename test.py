import numpy as np
import matplotlib.pyplot as plt

# Пример данных (замени на свои)
time = np.arange(0, 50, 0.01)  # Время от 0 до 50 с шагом 0.01
concentration = np.random.rand(len(time)) * np.exp(-((time - 25)**2) / (2 * 10**2))  # Пример концентрации

print(time, concentration, sep='\n')

# Связь времени с молекулярной массой (примем линейную зависимость)
k = 400  # Коэффициент, регулирующий масштаб масс (подбери по данным)
molecular_weights = time  # Переводим время в молекулярную массу

# Нормируем концентрацию
concentration /= np.sum(concentration)

# Построение графика ММР
plt.plot(molecular_weights, concentration)
plt.title('Результат расчета имитационной модели ММР')
plt.xlabel('Молекулярная масса')
plt.ylabel('Концентрация')
plt.show()