import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
 
#constants
a1 = 0.0134
b1 = 1
c1 = 4.35e-4
m1 = 1
a2 = 2.049
b2 = 0.563e-3
c2 = 0.528e5
m2 = 1
 
alpha0 = 0.05
alphaN = 0.01
l = 10
d = alphaN * l / (alphaN - alpha0)
cc = -alpha0 * d
 
T0 = 300
R = 0.5
Ft = 50
eps1 = 10e-4
eps2 = 10e-4
 
def plot2d(x, y, xlabel, ylabel, title):
    plt.plot(x, y, 'c')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.show()
 
def plot3d(x, y, z):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    #ax.set_title("")
    fig.show()
    plt.show()
 
#ax[n-1]+bx[n]+cx[n+1]=d
def tridiagonal_method(a, b, c, d):
    n = len(a) - 1;
    xi = [None for i in range(n+1)]
    xi[1] = -c[0]/b[0]
    eta = [None for i in range(n+1)]
    eta[1] = d[0]/b[0]
    for i in range(1, n):
        xi[i+1] = -c[i]/(a[i]*xi[i]+b[i])
        eta[i+1] = (d[i]-a[i]*eta[i])/(a[i]*xi[i]+b[i])
 
    res = [None for i in range(n+1)]
    res[n] =(d[n]-a[n]*eta[n])/(a[n]*xi[n]+b[n])
    for i in range(n-1, -1, -1):
        res[i] = xi[i+1]*res[i+1]+eta[i+1]
    return res
 
 
def k(T): 
    return a1 * (b1 + c1 * T**m1)
 
def c(T):
    return a2 + b2 * T**m2 - c2 / T**2
 
def alpha(x):
    return cc / (x - d)
 
def p(x):
    return 2/R * alpha(x)
 
#f(T) = f(x)
def f(x):
    return 2*T0/R * alpha(x)
 
 
 
def solve_equation_system(old_T, prev_T_by_time, tau, N):
    _a = []
    _b = []
    _d = []
    _f = []
    h = l / N
 
    def cal_plus_half(func, start, end):
        return (func(start) + func(end)) / 2
 
    def c_(n):
        return c(old_T[n])
 
    def c_plus(n):
        return cal_plus_half(c, old_T[n], old_T[n+1])
 
    def chi_plus(n):
        return cal_plus_half(k, old_T[n], old_T[n+1])
 
    def p_plus(n):
        return (p(n) + p(n+h)) / 2
 
    def f_plus(n):
        return (f(n) + f(n+h)) / 2
 
    def p_(n):
        return p(h*n)
 
    def f_(n):
        return f(h*n)
    
    #the left boundary condition x = 0
    _a.append(0)
    _b.append(h/8*c_plus(0)
              + h/4*c_(0)
              + tau/h*chi_plus(0)
              + tau*h/8*p_plus(0)
              + tau*h/4*p(0))
    _d.append(h/8*c_plus(0)
              - tau/h*chi_plus(0)
              + tau*h/8*p_plus(0))
    _f.append(h/8*c_plus(0)*(prev_T_by_time[0]+prev_T_by_time[1])
                + h/4*c_(0)*prev_T_by_time[0]
                + tau*Ft
                + tau*h/4*(f_plus(0) + f(0)))
 
    #1 <= n <= N-1
    for i in range(1, N):
        _a.append(tau/h*chi_plus(i-1))
        _d.append(tau/h*chi_plus(i))
        _b.append(-(_a[-1] + _d[-1] + c_(i)*h + p_(i)*h*tau))
        _f.append(-(f_(i)*h*tau + c_(i)*prev_T_by_time[i]*h))
 
    #the right boudary condition x = l
    _a.append(h/8*c_plus(N-1)
              - tau/h*chi_plus(N-1)
              + tau*h/8*p_plus(N-1))
    _b.append(h/4*c_(N)
              + h/8*c_plus(N-1)
              + tau*alphaN
              + tau/h*chi_plus(N-1)
              + tau*h/4*p_(N)
              + tau*h/8*p_plus(N-1))
    _d.append(0)
    _f.append(h/4*c_(N)*prev_T_by_time[N]
              + h/8*c_plus(N-1)*prev_T_by_time[N]
              + h/8*c_plus(N-1)*prev_T_by_time[N-1]
              + tau*alphaN*T0
              + tau*h/4*(f_(N)+f_plus(N-1)))
 
    return tridiagonal_method(_a, _b, _d, _f)
 
 
def check(res1, res2, eps):
    tmp = 0
    for i in range(len(res2)):
        tmp = max(tmp, math.fabs((res2[i]-res1[i])/res2[i]))
    if tmp <= eps:
        return True
    else:
        return False  
 
#simple iteration metho
def solve(tau, N):
    h = l / N
    time = 0
    old_T = [T0] * (N+1)
    new_T = [None] * (N+1)
    res = [old_T]
    while True:
        old_T = res[-1]
        while True:
            new_T = solve_equation_system(old_T, res[-1], tau, N)
            if check(old_T, new_T, eps1):
                break
            old_T = new_T
        time += tau
        res.append(new_T)
        if check(res[-2], res[-1], eps2):
              break
    return res
 
if __name__ == "__main__":
    tau = 1
    N = 100
    res = np.array(solve(tau, N))
 
    #3d surface T(x, t)
    x = np.linspace(0, l, N+1)
    y = np.arange(len(res))
    X, Y = np.meshgrid(x, y)
    Z = res
    plot3d(X, Y, Z)
 
    #Зависимость температуры от координаты стержня
    x = np.linspace(0, l, N+1)
    for T in res[:10]:
        plt.plot(x, T)
    plt.xlabel("Длина, см")
    plt.ylabel("Температура, K")
    plt.grid(True)
    plt.show()
 
    #Зависимость температуры от времени
    t = np.arange(len(res))
    for T in res.transpose()[:10]:
        plt.plot(t, T)
    plt.xlabel("Время, сек")
    plt.ylabel("Температура, K")
    plt.grid(True)
    plt.show()
