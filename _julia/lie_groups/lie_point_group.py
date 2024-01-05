from sympy import symbols, diff, Eq, Derivative

# Define the independent variables
x, y, t = symbols('x y t')

# Define the dependent variable and the differential equation
u = symbols('u', cls=Function)
eq = Eq(Derivative(u(x, y, t), t), u(x, y, t).diff(x, x) + u(x, y, t).diff(y, y))

# Define the infinitesimal generator, including the derivatives of the independent variables up to order 2
x1, y1, t1 = symbols('x1 y1 t1')
X = x1*Derivative(u, x) + y1*Derivative(u, y) + t1*Derivative(u, t) + symbols('X')

# Lie derivative
L_X = Derivative(X, u)

# Determine the determining equations
determining_eq = L_X(eq.lhs) - L_X(eq.rhs)

# Solve the determining equations
solutions = determining_eq.rhs.atoms(Symbol)

# The solutions of the determining equations will give you the form of the infinitesimal generator and the corresponding symmetries.
print(solutions)
