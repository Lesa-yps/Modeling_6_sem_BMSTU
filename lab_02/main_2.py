import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Константы
c = 3e10  # скорость света, см/с
R = 0.35  # радиус цилиндра, см
Tw = 2000  # температура стенки, K
T0 = 10000  # начальная температура, K
EPS = 1e-4  # точность метода стрельбы
P = 4

# Таблица значений k(T)
# Значения температуры и коэффициента поглощения, взятые из таблицы
# Они используются для интерполяции значений k(T)
t_values = np.array([2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000])
k_values = np.array([1.600E+00, 5.400E+00, 1.280E+01, 2.500E+01, 4.320E+01,
                      6.860E+01, 1.024E+02, 1.458E+02, 2.000E+02])

# Интерполяция k(T) в логарифмическом масштабе
# Используется логарифмическое преобразование для более точной интерполяции данных
log_k_interp = interp1d(np.log(t_values), np.log(k_values), kind='linear', fill_value='extrapolate')
def k_interp(T):
    """Вычисляет интерполированное значение k(T) с использованием логарифмической интерполяции."""
    return np.exp(log_k_interp(np.log(T)))

# Функция температуры
T = lambda r: T0 + (Tw - T0) * (r / R) ** P  # Линейное изменение температуры в цилиндре

# Функция Планка
up = lambda r: 3.084e-4 / (np.exp(4.799e4 / T(r)) - 1)  # Вычисление объемной плотности энергии

# Система ОДУ
def ode_system(r, u, F):
    """Формулирует систему дифференциальных уравнений для метода Рунге-Кутта."""
    k_r = k_interp(T(r))  # Вычисление коэффициента поглощения
    du_dr = - F * 3 * k_r / c #if r != 0 else F / (2 * k_r)  # Производная u по r (учет деления на ноль)
    if r == 0:
        dF_dr = c * k_r * (up(r) - u) / 2
    else:
        dF_dr = c * k_r * (up(r) - u) - F / r # Производная F по r
    return du_dr, dF_dr

# Метод Рунге-Кутта 4-го порядка
def runge_kutta(u0, F0, h=1e-2):
    """Численное решение системы ОДУ методом Рунге-Кутта 4-го порядка."""
    r_vals = np.arange(0, R + h, h)  # Создание массива значений r
    u_vals, F_vals = [u0], [F0]  # Начальные условия
    
    for r in r_vals[:-1]:
        u, F = u_vals[-1], F_vals[-1]  # Последние рассчитанные значения
        
        # Вычисление коэффициентов метода Рунге-Кутта
        k1, q1 = ode_system(r, u, F)
        k2, q2 = ode_system(r + h / 2, u + h * k1 / 2, F + h * q1 / 2)
        k3, q3 = ode_system(r + h / 2, u + h * k2 / 2, F + h * q2 / 2)
        k4, q4 = ode_system(r + h, u + h * k3, F + h * q3)
        
        # Итоговые значения на следующем шаге
        u_next = u + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        F_next = F + (h / 6) * (q1 + 2 * q2 + 2 * q3 + q4)
        
        u_vals.append(u_next)
        F_vals.append(F_next)
    
    return r_vals, np.array(u_vals), np.array(F_vals)

# Метод стрельбы (метод бисекции)
def shooting_method(chi_min=0.01, chi_max=1.0):
    """Подбор оптимального начального условия u(0) методом половинного деления."""
    while abs(chi_max - chi_min) > EPS:
        chi_mid = (chi_min + chi_max) / 2  # Среднее значение chi
        u0 = chi_mid * up(0)  # Начальное условие для u
        _, u_vals, F_vals = runge_kutta(u0, 0)  # Решение системы
        
        # Проверка условия 
        if (- F_vals[-1] + 0.39 * c * u_vals[-1]) > 0:
            chi_max = chi_mid
        else:
            chi_min = chi_mid
    
    return chi_mid

# Решение задачи
chi_opt = shooting_method()  # Поиск оптимального параметра chi (методом стрельбы)
u0_opt = chi_opt * up(0)  # Вычисление оптимального начального условия
r_vals, u_vals, F_vals = runge_kutta(u0_opt, 0)  # Численное решение системы методом Рунге-Кутта

# Построение графиков
plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.plot(r_vals, u_vals, 'b-', linewidth=2)
plt.title(f'Решение u(r) при xsi = {chi_opt:.6f}')
plt.xlabel('r')
plt.ylabel('u(r)')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(r_vals, F_vals, 'r-', linewidth=2)
plt.title('Функция F(r)')
plt.xlabel('r')
plt.ylabel('F(r)')
plt.grid(True)

plt.tight_layout()
plt.savefig('solution_2.png')

# Вывод оптимального значения chi
print(f'Оптимальное значение chi: {chi_opt}')
