export solve_y_xx


function solve_y_xx(;order=2, prnt = 0)
@variables x y y_x y_xx
#y_x = Dx(y); y_xx = Dx(Dx(y))

# define differential equation set to 0 of form \phi(x, y, yx, yxx) = 0
phi = y_xx

# define independent variables
independent_vars = [x, y, y_x, y_xx]

# create symbolic infintesimal generators
@variables eta(x,y) xi(x,y)


# define restrictions
restrictions = Dict([y_xx => 0])

# define free variables that can vary within the range of the restrictions
free_var = y_x

print("\n Attempting to find determining equations \n")
determining_eqns = find_determining_equations(phi, eta, xi, independent_vars, restrictions, free_var, order)

### TO DO ###
# now that we have determining equations, we need methods to find the actual values of the Infinitesimals
# given the determining equations
# use method 'two_var_power_series(expr=0, series_order,var1, var2)' to find power series of eta and xi
print("\n Attempting to solve using power series of infinitesimals \n")

xi_guess = 0; eta_guess = 0;
series_order = 3;
infin_gen = [xi, eta]
infin_gen_guess = [xi_guess, eta_guess]
variables = [x, y]

xi_sol, eta_sol = power_series_to_solve(infin_gen, infin_gen_guess, determining_eqns, series_order, variables)


return xi_sol, eta_sol
end