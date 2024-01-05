export two_var_power_series, wrapper_two_var_power_series, power_series_to_solve

function two_var_power_series(expr, series_order,var1, var2, const_term)
    
    c = 1
    for i in  0:series_order
        for j in 0:series_order
            if (i+j<=series_order)
            expr = expr + const_term[c]*var1^i * var2^j
            c = c + 1
            end
        end
    end

    return expr
end

function wrapper_two_var_power_series(expr1, expr2, series_order, var1, var2)
    
    nvar = 2 # set because non-general, two-variable power series
    nconst = binomial(series_order+nvar,nvar)
    @variables a[1:nconst] b[1:nconst]

    expr1 = two_var_power_series(expr1, series_order,var1, var2, a)
    expr2 = two_var_power_series(expr2, series_order,var1, var2, b)

    return expr1, expr2
end

function power_series_to_solve(infin_gen, infin_gen_guess, determining_eqns, series_order, variables)

    # develop power series for guess of xi and eta
xi_guess, eta_guess = wrapper_two_var_power_series(infin_gen_guess[1], infin_gen_guess[2], series_order, variables[1], variables[2])

# insert into determining equations
guesses = Dict([infin_gen[1] => xi_guess, infin_gen[2] => eta_guess])
free_var = [variables[1], variables[2]]
terms_to_zero = Array{Any}(undef, size(determining_eqns,1)) 

# collect like terms that need to satisfy equation, to send terms of dependent and independent variables to 0
for i in eachindex(determining_eqns)
    determining_eqns[i] = substitute(determining_eqns[i], guesses)
    #determining_eqns[i] = Symbolics.diff2term(Symbolics.value(determining_eqns[i]))
    determining_eqns[i] = expand_derivatives(determining_eqns[i]);

    determining_eqns[i] = simplify(determining_eqns[i]; expand=true,
                    threaded=false,
                    thread_subtree_cutoff=100,
                    rewriter=nothing)

    terms_to_zero[i] = collect_powers_two_var(determining_eqns[i], variables[1], variables[2], 0:series_order)
end


# seperate terms and send to zero by forming a dictionary, send those with more than one variable to each other, 
# unless it appears twice with different coefficients, then send to zero

sets_to_zero = Array{Any}(undef, size(terms_to_zero)) 
dict_to_merge = Dict([] => 0)
for l in range(0, 3)
    for i in 1:size(terms_to_zero, 1)
        for k in eachindex(terms_to_zero[i])

            terms_to_zero[i][k] = substitute(terms_to_zero[i][k], dict_to_merge)
            vars = Symbolics.get_variables(terms_to_zero[i][k])
            if size(vars,1) < 2
                for p in eachindex(vars)
                    sets_to_zero = Dict(vars[p]=> 0)
                        dict_to_merge = merge(dict_to_merge, sets_to_zero)
                end
            else # check if same factor for the two 
                term = copy(terms_to_zero[i][k])
                sets_to_zero = Dict(vars[p] => 1 for p in eachindex(vars))
                #print(vars)
                #print('\n')
                #print(term)
                #print('\n')

                term = substitute(term, sets_to_zero)
                term = substitute(term, dict_to_merge)
                #print(term)
                #print('\n')
                if term == 0    
                    sets_to_zero = Dict(vars[1] => vars[2]) #for p in eachindex(vars))
                    dict_to_merge = merge(dict_to_merge, sets_to_zero)
                elseif term != 0
                    sets_to_zero = Dict(vars[p] => 0 for p in eachindex(vars))
                    dict_to_merge = merge(dict_to_merge, sets_to_zero)
                end
            end
        end
    end
end



xi_guess = substitute(xi_guess, dict_to_merge)
eta_guess = substitute(eta_guess, dict_to_merge)

return xi_guess, eta_guess
end
