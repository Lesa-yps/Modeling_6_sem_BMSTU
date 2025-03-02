import numpy as np
import matplotlib.pyplot as plt

def Analytic_method(y):
    return -y**2/2 - 1/2 + np.exp(y**2)

def Picar_method(y):
    return 1/2 + y**2/2 + y**4/2 + y**6/6 + y**8/24 + y**10/240

# диапазон значений y
y_values = np.linspace(-1.7, 1.7, 1000)

# вычисление значения функций
picar_values = Picar_method(y_values)
analytic_values = Analytic_method(y_values)

# Построение графиков
plt.figure(figsize=(8, 6))
plt.plot(picar_values, y_values, label="Метод Пикара", color='blue')
plt.plot(analytic_values, y_values, label="Аналитическое решение", color='green')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()
plt.title("Сравнение метода Пикара и аналитического решения")
plt.show()