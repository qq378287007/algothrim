import numpy as np
from numpy import matrix as mat
from matplotlib import pyplot as plt

# 拟合函数
# A是系数，x是变量
def Func(Coefficient, x):
    return Coefficient[0] + (Coefficient[1] * x + Coefficient[2]) / (x * x + Coefficient[3] * x + Coefficient[4])  + \
        (Coefficient[5] * x + Coefficient[6]) / (x * x + Coefficient[7] * x + Coefficient[8])

# 函数偏导
# A是系数，x是变量，n是求偏导系数的序号
def Deriv(Coefficient, x, n):
    eps = 0.000001
    x1 = Coefficient.copy()
    x2 = Coefficient.copy()
    x1[n] -= eps
    x2[n] += eps
    y1 = Func(x1, x)
    y2 = Func(x2, x)
    dy = y2 - y1
    dx = eps * 2
    return dy / dx

data = np.loadtxt(open("ChinaSteel_35CS250H.csv", "rb"), delimiter=",") 
H = data[:, 0] #X
B = data[:, 1] #Y
# 数据个数
data_number = len(H)

# 系数个数
coefficient_number = 9

# 雅克比矩阵
J = mat(np.zeros((data_number, coefficient_number)))

# f(x)误差
fx = mat(np.zeros((data_number, 1))) 
fx_tmp = fx

# 系数初始化
Coefficient = mat(np.ones((coefficient_number, 1)))

conve = 100
u = 1
v = 2
lase_mse = 0
step = 0
while (conve):
    #计算f(x)和雅克比矩阵
    mse = 0
    for i in range(data_number):
        fx[i] =  Func(Coefficient, H[i]) - B[i] 
        mse += fx[i]**2
        for j in range(coefficient_number): 
            J[i, j] = Deriv(Coefficient, H[i], j) # 数值求导                                                    
    mse /= data_number  # 范围约束
    
    Hysen = J.T * J + u * np.eye(coefficient_number)
    dx = -Hysen.I * J.T * fx 
    xk_tmp = Coefficient.copy()
    xk_tmp += dx
    
    mse_tmp = 0
    for i in range(data_number):
        fx_tmp[i] =  Func(xk_tmp, H[i]) - B[i]  
        mse_tmp += fx_tmp[i, 0]**2
    mse_tmp /= data_number

    q = (mse - mse_tmp) / (0.5 * dx.T * (u * dx - J.T * fx)) 
    q = q[0, 0]
    if q > 0:
        mse = mse_tmp
        max_value = max(1.0/3.0, 1 - pow(2 * q - 1, 3))
        u *= max_value
        v = 2
    else:
        u *= v
        v *= 2
    Coefficient = xk_tmp
    
    error = abs(mse - lase_mse)
    step += 1
    print("step = %d, abs(mse - lase_mse) = %.8f" %(step, error))
    if error < 0.000001:
        break

    lase_mse = mse
    conve -= 1
print(Coefficient)

plt.figure(0)
plt.scatter(H, B, s = 4)

a = np.zeros(coefficient_number)
for i in range(coefficient_number):
    a[i] = Coefficient[i, 0]
print(a)

Z = [Func(a, x) for x in H] #用拟合好的参数画图
plt.plot(H, Z, 'r')
plt.show()
