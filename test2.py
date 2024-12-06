# Пример массива
numbers = [3, 15, 8, 22, 5]

# Шаги
steps = 100

# Номер шага (например, 25-й шаг)
step_number = 25

# Вычисление
min_value = min(numbers)
max_value = max(numbers)
step = (max_value - min_value) / steps

# Число, соответствующее номеру шага
result = min_value + (step_number * step)

print((15 - min_value) / step)