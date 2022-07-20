from math import ceil, sqrt

def output(s):
    if type(s) == float:
        if s > 1000000:
            return '{:.8e}'.format(s)
        return '{:.8f}'.format(s)
    elif type(s) == int:
        return str(s)
    else:
        return s

def f(x, u):
    return x ** 2  + u ** 2

# Эйлера
def euler(n, h, x, y):
    y_out = []
    for i in range(n):
        try:
            y += h * f(x, y)
            y_out.append(y)
            x += h
        except OverflowError:
            y_out.append('overflow')
            for j in range(i, n-1):
                y_out.append('-----')
            break
    return y_out

# РунгеКутта
def rungeKutta(n, h, x, y):
    y_out = [y]
    alpha = 1 / 2
    h1 = h / (2 * alpha)
    for i in range(n):
        try:
            k1 = f(x, y)
            k2 = f(x + h1, y + h1 * k1)
            y += h * ((1 - alpha) * k1 + alpha * k2)
            y_out.append(y)
            x += h
        except:
            y_out.append('overflow')
            for j in range(i, n-1):
                y_out.append('-----')
            break
    return y_out
        

# Пикар
def picar(n, h, x, y0):
    def f1(a):
        return a ** 3 / 3
    def f2(a):
        return f1(a) + a ** 7 / 63
    def f3(a):
        return f2(a) +  (a ** 11) * (2 / 2079) + (a ** 15) / 59535
    def f4(a):
        return f3(a) + (a ** 15)*(2 / 93555) + (a ** 19)*(2 / 3393495) + (a ** 19)*(2 / 2488563) + \
    (a ** 23)*(2 / 86266215) + (a ** 23)*(1 / 99411543) + (a ** 27)*(2 / 3341878155) + (a ** 31)*(1 / 109876902975)

    y_out = [[y0, y0, y0, y0]]
    for i in range(n-1):
        x += h
        y_f1 = f1(x)
        y_f2 = f2(x)
        y_f3 = f3(x)
        y_f4 = f4(x)
        y_out.append([y_f1, y_f2, y_f3, y_f4])
    return y_out
        

def main():
    h = 10 ** -5 
    x = 0
    y0 = 0
    end = 2.1

    n = ceil(abs(end - x)/h)+1

    x_arr = [x + h*i for i in range(n)]
    y_picar = picar(n, h, x, y0)
    y_euler = euler(n, h, x, y0)
    y_RungeKutta = rungeKutta(n, h, x, y0)

    print("|    x    |    Пикара 1   |    Пикара 2   |    Пикара 3   |    Пикара 4   |     Эйлер     |    РунгеКутт |")
    print("-"*75)
    output_step = int(n/h) # выводим только 100 значений в таблице 
    for i in range(int(1.98 / h), n, 1):
        print("|{:^9.5f}|{:^15.8f}|{:^15.8f}|{:^15.8f}|{:^15.8f}|{:^15s}|{:^15s}|".format(x_arr[i], y_picar[i][0], y_picar[i][1], y_picar[i][1], y_picar[i][1], output(y_euler[i]), output(y_RungeKutta[i])))


main()
