import sympy as sp
from tabulate import tabulate
import matplotlib.pyplot as plt

import sys
sys.stdout.reconfigure(encoding='utf-8')

# константы
H = 0.001  # шаг
X_MAX = 1.64 # 1.64724785  # граница расчётов по x

# МЕТОД ПИКАРА

# определение символов
x = sp.symbols('x')
u = sp.Function('u')(x)

# явное вычисление приближений Пикара (1-4)
u1 = sp.integrate(x + 0**3, (x, 0, x))
u2 = sp.integrate(x + (u1)**3, (x, 0, x))
u3 = sp.integrate(x + (u2)**3, (x, 0, x))
u4 = sp.integrate(x + (u3)**3, (x, 0, x))

# преобразование в численные функции
Picar_u1_func = sp.lambdify(x, u1, 'numpy')
Picar_u2_func = sp.lambdify(x, u2, 'numpy')
Picar_u3_func = sp.lambdify(x, u3.simplify(), 'numpy')
Picar_u4_func = sp.lambdify(x, u4.simplify(), 'numpy')

# МЕТОД ЭЙЛЕРА
def Euler_func(x, y, h, n):
    def f(x, y):
        return x + y**3
    arr_res = list()
    arr_res.append((x, y))
    for _ in range(n):
        y += h * f(x, y)
        x += h
        arr_res.append((x, y))
    return arr_res

# создание таблицы
n = int(X_MAX / H)
x_values = [i * H for i in range(n + 1)]
picard_values = [[x, Picar_u1_func(x), Picar_u2_func(x), Picar_u3_func(x), Picar_u4_func(x)] for x in x_values]
euler_values = Euler_func(0, 0, H, n)

# объединение данных в одну таблицу
combined_values = []
for i in range(len(x_values)):
    combined_values.append([x_values[i], *picard_values[i][1:], euler_values[i][1]])

table = tabulate(combined_values, headers=["x", "Picard 1", "Picard 2", "Picard 3", "Picard 4", "Euler"], tablefmt="grid")

print("Table of values ​​of Picard and Euler methods:")
print(table)


# построение графика
plt.figure(figsize=(10, 6))
plt.plot(x_values, [picard_value[1] for picard_value in picard_values], label='Picard 1')
plt.plot(x_values, [picard_value[2] for picard_value in picard_values], label='Picard 2')
plt.plot(x_values, [picard_value[3] for picard_value in picard_values], label='Picard 3')
plt.plot(x_values, [picard_value[4] for picard_value in picard_values], label='Picard 4')
plt.plot(*zip(*euler_values), label='Euler')

plt.xlabel("x")
plt.ylabel("u(x)")
plt.title("Comparison of Picard and Euler Methods")
plt.legend()
plt.grid()
plt.savefig("3_graph.png", dpi=300) 
plt.show()
