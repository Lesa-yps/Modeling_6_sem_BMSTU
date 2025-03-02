H = 0.5  # начальный шаг
CONTROL_POINT = 1.5  # длина x
EPS = 10**(-7)  # точность

def f(x, y):
    return x + y**3

def calc_n(x1, x2, h):
    return int((x2 - x1) / h) + 1

def Euler_func(x, y, h):
    while 1:
        try:
            y += h * f(x, y)
        except Exception:
            print("Too large x =", x, "y =", y, "h =", h)
            break
        x += h
    return x, y


def Euler_func_xmax(x, y, h, xmax):
    while x < xmax:
        try:
            y += h * f(x, y)
        except Exception:
            print("Error: too large x =", x - h, "y =", y, "h =", h)
            break
        x += h
    return y

def adaptive_Euler(x0, y0, h, control_point, eps):
    y_h2 = Euler_func_xmax(x0, y0, h, control_point)
    
    while True:
        y_h = y_h2
        h /= 2
        print("STEP =", h)
        y_h2 = Euler_func_xmax(x0, y0, h, control_point)
        
        if abs((y_h - y_h2) / y_h2) < eps:
            print(f"FINAL STEP = {h}\n")
            break

    return h

h = adaptive_Euler(0, 0, H, CONTROL_POINT, EPS)
x_max, y_max = Euler_func(0, 0, h)
print(f"x_max: {x_max}, y_max: {y_max}, step: {h}")