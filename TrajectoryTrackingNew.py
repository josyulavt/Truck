import pydrake
import pydrake.symbolic as sym
from pydrake.all import (Variable, OsqpSolver, Solve, SnoptSolver, eq, IpoptSolver)
from pydrake.solvers import MathematicalProgram 
import numpy as np

def continuous_kinematics(x, u):
    # x = [x, y, theta0, theta1, v, delta]
    # u = [v, delta]
    d1 = 1
    m = np
    x1dot = v * m.cos(theta0 - theta1) * m.cos(theta1)
    y1dot = v * m.cos(theta0 - theta1) * m.sin(theta1)
    theta0dot = (v / d1) * m.tan(delta)
    theta1dot = (v / d1) * m.sin(theta0 - theta1)
    x_d = np.array([
        x1dot,
        y1dot,
        theta0dot,
        theta1dot,
        u[0],
        u[1]
    ], dtype=np.float64)
    return x_d

def lqr(x_bar, u_bar, dt):
    d1 = 1
    # initialize problem
    Q = np.eye(5, dtype=np.float64)
    R = np.eye(2, dtype=np.float64)

    prog = MathematicalProgram()
    x = prog.NewContinuousVariables(6, 2, "x")
    u = prog.NewContinuousVariables(2, "u")

    A = np.array([
        [1, 0, 0, 0, np.cos(x_bar[2] - x_bar[3]), -u_bar[0] * np.sin(x_bar[2] - x_bar[3])],
        [0, 1, 0, 0, np.sin(x_bar[2] - x_bar[3]), u_bar[0] * np.cos(x_bar[2] - x_bar[3])],
        [0, 0, 1, 0, 0, u_bar[0] / (d1 * np.cos(u_bar[1]) ** 2)],
        [0, 0, 0, 1, -np.cos(x_bar[2] - x_bar[3]), -u_bar[0] * np.sin(x_bar[2] - x_bar[3])],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]
    ], dtype=np.float64)

    B = np.array([
        [0, 0],
        [0, 0],
        [0, 0],
        [0, 0],
        [np.cos(x_bar[2] - x_bar[3]) / d1, 0],
        [0, 1 / d1]
    ], dtype=np.float64)

    A_x = A.dot(x_bar)
    B_u = B.dot(u_bar)

    dx = x[:, 0] - x_bar
    du = u[:] - u_bar
    dxdt = A_x.dot(dx) + B_u.dot(du)
    fxu = continuous_kinematics(x[:, 0], u)
    for i in range(6):
        prog.AddConstraint(x[i, 1] == x[i, 0] + (fxu[i] + dxdt[i]) * dt)
    prog.AddQuadraticErrorCost(Q, x_bar, x[:, 0])
    prog.AddQuadraticErrorCost(R, u_bar, u)

    solver = IpoptSolver()

    result = solver.Solve(prog)
    x_opt = result.GetSolution(x)
    return x_opt

# test lqr with fake data

lqr(np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1]), np.array([0, 0]), 0.1)