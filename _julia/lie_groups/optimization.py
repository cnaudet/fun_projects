import pyomo.environ as pyo

m = pyo.ConcreteModel()

d = 20

m.PRODUCTS = pyo.Set(initialize=["U","V"])
m.x = pyo.Var(domain=pyo.NonNegativeReals)
m.y = pyo.Var(m.PRODUCTS, domain=pyo.NonNegativeReals)


@m.Constraint()
def demand_ub(m):
    return m.y["U"] <= 40

@m.Constraint()
def demand_lb(m):
    return m.y["U"]>=d

@m.Constraint()
def labor_A(m):
    return m.y["U"] + m.y["V"] <= 80

@m.Constraint()
def labor_B(m):
    return 2*m.y["U"] + m.y["V"] <= 100

@m.Constraint()
def production(m):
    return m.x >= 10*m.y["U"] + 9*m.y["V"]

revenue = 279*m.y["U"] + 210*m.y["V"]
expenses = 130*m.y["U"] + 90*m.y["V"] + 10*m.x

@m.Objective(sense=pyo.maximize)
def profit(m): 
    return revenue - expenses


solver = pyo.SolverFactory('glpk')
solver.solve(m).write()
