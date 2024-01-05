
export find_determining_equations, construct_group_operator, construct_infinitesimals, total_differentiation_op, collect_powers

#struct Differential <: Function
#end


function find_determining_equations(phi, eta, xi, independent_vars, restrictions, free_var, order)

    # construct higher order infinitesimals up to order `order'
    # return infinitesimal derivatives infin_der = [xi, eta, eta_1, eta_2]
    print("\n Constructing Infinitesimals \n")
    infin_der = construct_infinitesimals(eta, xi, independent_vars, order)

    # simplify
    infin_der = simplify(infin_der; expand=true,
                       threaded=false,
                       thread_subtree_cutoff=100,
                       rewriter=nothing)

    # redefine because it refines above
    #independent_vars = [x, y, y_x, y_xx]
    # construct group operator X applied to original PDE phi, as passed into first argument

    print("\n Constructing Group Operators \n")
    X = construct_group_operator(phi, independent_vars, infin_der)

    # begin solving for infinitesimal generators by substituting and setting terms = 0
    # simplify expression
    X = expand_derivatives(X);

    # substitute restrictions
    X = substitute(X, restrictions)

    X = simplify(X; expand=true,
                      threaded=false,
                       thread_subtree_cutoff=100,
                    rewriter=nothing)

    X_latex = latexify(X)

    determining_eqns = collect_powers(X, free_var, 0:3)
return determining_eqns
end

function construct_group_operator(phi, independent_vars, infin_der)

    diffs  = Vector{Any}(undef, size(independent_vars,1)) 
    for i in eachindex(independent_vars)
        diffs[i] = Differential(independent_vars[i])
    end

    X = 0
    for i in 1:size(independent_vars,1)
        X = X + infin_der[i]diffs[i](phi)
    end
    return X
end

function construct_infinitesimals(eta, xi, independent_vars, p::Integer)

# independent variables of form x, y, y_x, y_xx, y_xxx, ...

    Dx = Differential(independent_vars[1])
    infin_der = Vector{Any}(undef, p+1) 
    infin_der[1] = eta[1]
    for i in 2:p+1
        infin_der[i] = total_differentiation_op(infin_der[i-1], independent_vars) - independent_vars[i+1]total_differentiation_op(xi, independent_vars)

    end
    infin_der = [xi; infin_der]
    return infin_der
end

function total_differentiation_op(F, independent_vars)
    # total differentiation operator with respect to the jth independent variable
    # F is function to differentiate with respect to
    # independent variables of form x, y, y_x, y_xx, y_xxx, ...

    # form derivitives with respect to independent variables, term 1, 2, ... p
    Dx = Differential(independent_vars[1])
    diffs  = Vector{Any}(undef, size(independent_vars,1)) 
    p = size(independent_vars, 1)
    for i in 1:size(independent_vars,1)
        diffs[i] = Differential(independent_vars[i])
    end
    
    final = Dx(last(independent_vars))
    scaling_fact = Vector{Any}(undef, p) 
    scaling_fact[1] = 1
    scaling_fact[2:p] = [independent_vars[3:p]; final]

    D  = 0
    for k in 1:size(scaling_fact,1)
        D = D + scaling_fact[k]diffs[k](F)
    end

    return D
end

function collect_powers(eq, x_p, ns; max_power=100)
    eq = substitute(expand(eq), Dict(x_p^j => 0 for j=last(ns)+1:max_power))

    eqs = []
    zero_order_term = 0
    for i in ns
        powers = Dict(x_p^j => (i==j ? 1 : 0) for j=1:last(ns))
        push!(eqs, substitute(eq, powers))
        (i==0 ? zero_order_term = eqs[1] : eqs[i+1] = eqs[i+1] - zero_order_term)
    end
    eqs
end

function collect_powers_two_var(eq, x_p, y_p, ns; max_power=100)
    eq = substitute(expand(eq), Dict(x_p^j => 0 for j=last(ns)+1:max_power))

    eqs = []
    zero_order_term = 0
    k_order0_term = 0
    i_order0_term = zeros(last(ns)+1)*x_p
    c = 0
    for k in ns
        for i in ns
            
            powers3 = Dict((x_p^j)y_p^k => (i==j ? 1 : 0) for j=1:last(ns))
            powers1 = Dict((x_p^j) => (i==j ? 1 : 0) for j=1:last(ns))
            powers2 = Dict(y_p^j => (k==j ? 1 : 0) for j=1:last(ns))
            subs = merge(powers1, powers2, powers3)
            push!(eqs, substitute(eq, subs))
#            push!(eqs, substitute(eq, powers1))
            (i==0&&k==0 ? zero_order_term = eqs[c+1] : eqs[c+1] = eqs[c+1] - zero_order_term)
            (k==0&&i>0 ? i_order0_term[i+1] = eqs[c+1] : eqs[c+1] = eqs[c+1] - i_order0_term[i+1])
            (i==0&&k>0 ? k_order0_term = eqs[c+1] : eqs[c+1] = eqs[c+1] - k_order0_term)
            c = c + 1
        end
    end
    eqs
end
 





