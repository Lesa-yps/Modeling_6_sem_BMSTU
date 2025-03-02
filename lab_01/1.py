import matplotlib.pyplot as plt

STEP = 1e-3
DIST = 100

def Euler_func(x, y, z):
    return (-2 * z - y) / (4 * x)

def Euler_method(x, y, z, h, n):
    arr_res = list()
    arr_res.append((x, y))
    for _ in range(n):
        y += h * z
        z += h * Euler_func(x, y, z)
        x += h
        arr_res.append((x, y))
    return arr_res

def Taylor_func(x):
    return 1 - x/2 + x**2/24 - x**3/720 + x**4/40320

def Taylor_series(x, h, n):
    y = Taylor_func(x)
    arr_res = list()
    arr_res.append((x, y))
    for _ in range(n):
        x += h
        y = Taylor_func(x)
        arr_res.append((x, y))
    return arr_res

if __name__ == "__main__":
    x = 10**(-3)
    y = 1
    z = -1/2
    h = STEP
    n = int(DIST / h)
    
    # получение результатов вычислений
    arr_Euler = Euler_method(x, y, z, h, n)
    arr_Taylor = Taylor_series(x, h, int(n / 20))

    # разделение x и y для удобства построения графиков
    x_vals_Euler, y_vals_Euler = zip(*arr_Euler)
    x_vals_Taylor, y_vals_Taylor = zip(*arr_Taylor)

    # отрисовка графиков
    plt.figure(figsize=(8, 5))
    plt.plot(x_vals_Euler, y_vals_Euler, label='Метод Эйлера', marker='o', markersize=2, linestyle='-', alpha=0.7)
    plt.plot(x_vals_Taylor, y_vals_Taylor, label='Ряд Тейлора', marker='x', markersize=2, linestyle='-', alpha=0.7)

    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('Сравнение метода Эйлера и ряда Тейлора')
    plt.legend()
    plt.grid()
    plt.show()
