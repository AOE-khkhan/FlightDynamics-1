import numpy as np

def f(t, x, **params):

    #a = params['a']
    #c = params['c']

    #f1 = a * (x[0] - x[0] * x[1])
    #f2 = -c * (x[1] - x[0] * x[1])
    f1 = xdot[0]
    f2 = xdot[1]
    f3 = xdot[2]
    f4 = xdot[3]
    f5 = xdot[4]
    f6 = xdot[5]
    f7 = xdot[6]
    f8 = xdot[7]
    f9 = xdot[8]
    f10 = xdot[9]
    f11 = xdot[10]
    f12 = xdot[11]


    return np.array([f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12], dtype=np.float)

def df(t, x, **params):

    eps = 1e-10
    J = np.zeros([len(x), len(x)], dtype=np.float)

    for i in range(len(x)):
        x1 = x.copy()
        x2 = x.copy()

        x1[i] += eps
        x2[i] -= eps

        f1 = f(t, x1, **params)
        f2 = f(t, x2, **params)

        J[:, i] = (f1 - f2) / (2 * eps)

    return J

t = 0
xdot = [-9.30762134e-11,0.00000000e+00,3.19822391e-10,9.00285110e+01,0.00000000e+00,-3.29685917e-16,0.00000000e+00,-9.78779587e-11,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00]
x= np.array(list(range(1,13)), dtype=np.float)
print(df(t, x))

#t = 0
#x = np.array([1,2], dtype=np.float)
#print(df(t, x, a=1, c=1))



# class ODEModel:

#     def __init__(self, eps=1e-10):
#         self.eps = eps
#         self.farray  = []
    
#     def add_function(self, f):
#         self.farray.append(f)
    
#     def f(self, t, x):
#         return np.array([fi(t,x) for fi in self.farray], dtype=np.float)
    
#     def df(self, t, x):
#         J = np.zeros([len(x), len(x)], dtype=np.float)

#         for i in range(len(x)):
#             x1 = x.copy()
#             x2 = x.copy()

#             x1[i] += self.eps
#             x2[i] -= self.eps

#             f1 = self.f(t, x1)
#             f2 = self.f(t, x2)

#             J[:, i] = (f1 - f2) / (2 * self.eps)

#         return J


# F = ODEModel(eps = 1e-12)

# # eq1 = lambda t,u : u[1]
# # eq2 = lambda t,u : u[2]
# # eq3 = lambda t,u : u[3]
# # eq4 = lambda t,u : -8 * u[0] + np.sin(t) * u[1] - 3 * u[3] + t**2

# # F.add_function(eq1)
# # F.add_function(eq2)
# # F.add_function(eq3)
# # F.add_function(eq4)

# # t = 0
# # x = np.array([1,2,3,4], dtype=np.float)
# # print(F.df(t,x))

# xdot = [-9.30762134e-11,0.00000000e+00,3.19822391e-10,9.00285110e+01,0.00000000e+00,-3.29685917e-16,0.00000000e+00,-9.78779587e-11,0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00]

# eq1 = xdot[0]
# eq2 = xdot[1]
# eq3 = xdot[2]
# eq4 = xdot[3]
# eq5 = xdot[4]
# eq6 = xdot[5]
# eq7 = xdot[6]
# eq8 = xdot[7]
# eq9 = xdot[8]
# eq10 = xdot[8]
# eq11 = xdot[9]
# eq12 = xdot[10]


# F.add_function(eq1)
# F.add_function(eq2)
# F.add_function(eq3)
# F.add_function(eq4)
# F.add_function(eq5)
# F.add_function(eq6)
# F.add_function(eq7)
# F.add_function(eq8)
# F.add_function(eq9)
# F.add_function(eq10)
# F.add_function(eq11)
# F.add_function(eq12)

# t = 0
# x = np.array([1,2,3,4,5,6,7,8,9,10,11], dtype=np.float)
# print(F.df(t,x))
