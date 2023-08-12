import pydrake
import pydrake.symbolic as sym
from pydrake.all import (Variable, OsqpSolver, Solve, SnoptSolver, eq, IpoptSolver)
from pydrake.solvers import MathematicalProgram 
import numpy as np

def continuous_kinematics(x, u):
    #x = [x, y, theta0, theta1, v, delta]
    #u = [v, delta]
    d1 = 1
    print(x[0])
    # x = x[0]
    # y = x[1]
    theta0 = x[2]
    theta1 = x[3]
    v = x[4]
    delta = x[5]
    m = sym if x.dtype == object else np 
    x1dot = v*m.cos(theta0-theta1)*m.cos(theta1);
    y1dot = v*m.cos(theta0-theta1)*m.sin(theta1);
    theta0dot = (v/d1)*m.tan(delta); 
    theta1dot = (v/d1)*m.sin(theta0-theta1); 
    x_d = np.array([
        x1dot, 
        y1dot,
        theta0dot,
        theta1dot,    
        u[0],   
        u[1]       
    ])
    return x_d


def discrete_dynamics(x, u, dt):
    # forward Euler
    x_next = x + dt * continuous_kinematics(x, u)
    return x_next

# def linear_dynamics_constraints(prog, decision_vars, x_bar, u_bar, N, dt):
#     '''
#     Dynamics Constraints -- error state form
#     x_bar - nominal point; dx = (x - x_bar)
#     A, B, C = get_linear_dynamics(derivs, x_bar[:,n], u_bar[:,n])
#     '''
#     x, u, _, _ = decision_vars
#     derivs = dynamics.derivatives(dynamics.discrete_dynamics, n_x, n_u)
#     for n in range(N-1):
#         A, B = dynamics.get_linear_dynamics(derivs, x_bar[:, n], u_bar[:, n])
#         dx = x[:, n] - x_bar[:, n]
#         du = u[:, n] - u_bar[:, n]
#         dxdt = A@dx + B@du
#         fxu_bar = dynamics.car_continuous_dynamics(x_bar[:, n], u_bar[:, n])
#         for i in range(n_x):
#             prog.AddConstraint(x[i, n+1] == x[i, n] +
#                                (fxu_bar[i] + dxdt[i])*dt)
#     return prog

def lqr(x_bar, u_bar, dt): 
#intialize problem 
    Q = np.eye(6, dtype=np.float64)
    R = np.eye(2, dtype=np.float64)
    
    prog = MathematicalProgram()
    x = prog.NewContinuousVariables(6, 2, "x")
    u = prog.NewContinuousVariables(2, "u")

    x_sym = np.array([sym.Variable("x" + str(i)) for i in range(6)])
    u_sym = np.array([sym.Variable("u" + str(i)) for i in range(2)])
    
    A = sym.Jacobian(continuous_kinematics(x_sym, u_sym), x_sym)
    B = sym.Jacobian(continuous_kinematics(x_sym, u_sym), u_sym)
    
    A_x = sym.Evaluate(A, dict(zip(x_sym, x_bar)))
    B_u = sym.Evaluate(B, dict(zip(x_sym, x_bar)))
    
    dx = x[:,0] - x_bar
    du = u[:] - u_bar
    dxdt = A_x@dx + B_u@du
    # print(dxdt.shape)
    fxu = continuous_kinematics(x[:,0], u)
    for i in range(6): 
        prog.AddConstraint(x[i, 1] == x[i, 0] + (fxu[i] + dxdt[i])*dt)
    prog.AddQuadraticErrorCost(Q, x_bar, x[:,0])
    prog.AddQuadraticErrorCost(R, u_bar, u)
    
    solver = IpoptSolver()
    
    result = solver.Solve(prog)
    u_sol = result.GetSolution(u)
    x_sol = result.GetSolution(x)
    print("u_sol", u_sol)
    print("x_sol", x_sol)
    return u_sol, x_sol


#test lqr with fake data

# lqr(np.array([0.1,0.1,0.1,0.1,0.1,0.1]), np.array([0,0]), 0.1)
import matplotlib.pyplot as plt

x = [0, 0.0516223, 0.0891949, 0.0881861, 0.008462, -0.256918, -1.03396, -3.14422, -8.2997, -19.5812, -41.9364] 
y = [5, 5.08979, 5.16027, 5.124, 4.91303, 4.54727, 4.11682, 3.69573, 3.25676, 2.65966, 1.73982]
theta0 = [0.0] 
theta1 = [0.0]
v = [0.0] 
delta = [0.0]

x_act = [] 
y_act = []
v_act = []
delta_act = []

for i in range(10):
    x_bar = np.array([x[i], y[i], 0, 0, 0, 0])
    u_bar = np.array([0, 0])
    u_sol, x_sol = lqr(x_bar, u_bar, 0.1)
    x_new = discrete_dynamics(x_bar, u_sol, 0.1)
    x_act.append(x_new[0])
    y_act.append(x_new[1])
    theta0.append(x_new[2])
    theta1.append(x_new[3])
    v.append(x_new[4])
    delta.append(x_new[5])
    v_act.append(u_sol[0])
    delta_act.append(u_sol[1])

plt.plot(x_act, y_act, 'r-')    
plt.plot(x, y, 'b+')
plt.show()
 #fixed for tracking with lqr, need to implement tracking with mpc