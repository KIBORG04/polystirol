import numpy as np
import plotly.graph_objects as go

# Данные для графика
x = np.linspace(-10, 10, 100)
y = np.linspace(-10, 10, 100)
x, y = np.meshgrid(x, y)
z = np.sin(np.sqrt(x**2 + y**2))  # Функция для поверхности

print(z, len(x), len(y))
# Создание графика
fig = go.Figure(data=[go.Surface(z=z, x=x, y=y)])

# Настройка осей
fig.update_layout(scene=dict(
    xaxis_title="X Axis",
    yaxis_title="Y Axis",
    zaxis_title="Z Axis"
))

# Отображение
#fig.show()
