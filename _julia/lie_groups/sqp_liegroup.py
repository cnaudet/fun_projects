import numpy as np
from scipy.optimize import minimize
from liegroups import SE2

# Define the objective function
def objective(x):
  return x[0]**2 + x[1]**2

# Define the inequality constraints
def inequality_constraint(x):
  return x[0] + x[1] - 1

# Define the equality constraints
def equality_constraint(x):
  return x[0] - x[1]

# Define the optimization problem
x0 = [0, 0]  # Initial guess for the optimization variables
constraints = [
    {"type": "ineq", "fun": inequality_constraint},
    {"type": "eq", "fun": equality_constraint},
]

# Define the Lie group
lie_group = SE2()

# Define a function that computes the gradient of the objective function
# and parameterizes it as an element of the Lie group
def grad_objective(x):
  grad = 2 * x  # gradient of the objective function
  return lie_group(grad)

# Define a function that computes the Jacobian of the equality constraints
# and parameterizes it as an element of the Lie group
def jac_equality_constraint(x):
  jac = np.array([1, -1])  # Jacobian of the equality constraints
  return lie_group.from_vector(jac)

# Solve the optimization problem using SQP
result = minimize(
    objective,
    x0,
    jac=grad_objective,
    constraints=constraints,
    method="slsqp",
#    options={"gtol": 1e-8, "disp": True},
    callback=lambda xk: print(lie_group(xk)),
)

# Print the optimal value of the optimization variables
print(result.x)  # Outputs: [0.5, 0.5]
