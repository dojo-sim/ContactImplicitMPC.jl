const MOI = MathOptInterface

# user should define obj, con!
function moi_obj(x)
    @error "objective not defined"
end

function moi_con!(c, x)
    @error "constraints not defined"
end

function ∇moi_obj!(g, x)
    @error "objective gradient not defined"
end

function ∇moi_con!(j, x)
    @error "constraints Jacobian not defined"
end

# user can overwrite primal_bounds, constraint_bounds
function primal_bounds(n)
    x_l = -Inf * ones(n)
    x_u = Inf * ones(n)
    return x_l, x_u
end

function constraint_bounds(m; idx_ineq=(1:0))
    c_l = zeros(m)

    c_u = zeros(m)
    c_u[idx_ineq] .= Inf

    return c_l, c_u
end

struct ProblemMOI <: MOI.AbstractNLPEvaluator
    n_nlp::Int
    m_nlp::Int
    idx_ineq
    primal_bounds
    constraint_bounds
    sparsity_jac
    sparsity_hess
    hessian_lagrangian::Bool
end

function ProblemMOI(n_nlp, m_nlp;
        idx_ineq=(1:0),
        sparsity_jac=sparsity_jacobian(n_nlp, m_nlp),
        sparsity_hess=sparsity_hessian(n_nlp, m_nlp),
        primal_bounds=primal_bounds(n_nlp),
        constraint_bounds=constraint_bounds(m_nlp, idx_ineq=idx_ineq),
        hessian_lagrangian=false)

    ProblemMOI(n_nlp,m_nlp,
        idx_ineq,
        primal_bounds,
        constraint_bounds,
        sparsity_jac,
        sparsity_hess,
        hessian_lagrangian)
end

function MOI.eval_objective(prob::MOI.AbstractNLPEvaluator, x)
    moi_obj(x)
end

function MOI.eval_objective_gradient(prob::MOI.AbstractNLPEvaluator, grad_f, x)
    ∇moi_obj!(grad_f, x)
    return nothing
end

function MOI.eval_constraint(prob::MOI.AbstractNLPEvaluator, g, x)
    moi_con!(g, x)
    return nothing
end

function MOI.eval_constraint_jacobian(prob::MOI.AbstractNLPEvaluator, jac, x)
    ∇moi_con!(jac, x)
    return nothing
end

function row_col!(row,col,r,c)
    for cc in c
        for rr in r
            push!(row,convert(Int,rr))
            push!(col,convert(Int,cc))
        end
    end
    return row, col
end

function row_col_cartesian!(row,col,r,c)
    for i = 1:length(r)
        push!(row,convert(Int,r[i]))
        push!(col,convert(Int,c[i]))
    end
    return row, col
end

# user can overwrite sparsity_jacobian and sparsity_hessian
function sparsity_jacobian(n, m)

    row = []
    col = []

    r = 1:m
    c = 1:n

    row_col!(row,col,r,c)

    return collect(zip(row,col))
end

function sparsity_hessian(n,m)

    row = []
    col = []

    r = 1:m
    c = 1:n

    row_col!(row,col,r,c)

    return collect(zip(row,col))
end

function MOI.features_available(prob::MOI.AbstractNLPEvaluator)
    if prob.hessian_lagrangian
        return [:Grad, :Jac, :Hess]
    else
        return [:Grad, :Jac]
    end
end

MOI.initialize(prob::MOI.AbstractNLPEvaluator, features) = nothing
MOI.jacobian_structure(prob::MOI.AbstractNLPEvaluator) = prob.sparsity_jac
MOI.hessian_lagrangian_structure(prob::MOI.AbstractNLPEvaluator) = prob.sparsity_hess

function MOI.eval_hessian_lagrangian(prob::MOI.AbstractNLPEvaluator, H, x, σ, λ)
   @warn "Hessian of Lagrangian not implemented"
    return nothing
end

function solve(x0, prob::MOI.AbstractNLPEvaluator;
        nlp=:ipopt,
        tol=1.0e-5,
        c_tol=1.0e-5,
        max_iter=1000)

    x_l, x_u = prob.primal_bounds
    c_l, c_u = prob.constraint_bounds

    nlp_bounds = MOI.NLPBoundsPair.(c_l,c_u)
    block_data = MOI.NLPBlockData(nlp_bounds, prob, true)

    if nlp == :ipopt
        solver = Ipopt.Optimizer()
        solver.options["max_iter"] = max_iter
        solver.options["tol"] = tol
        solver.options["compl_inf_tol"] = tol
        solver.options["dual_inf_tol"] = tol
        solver.options["constr_viol_tol"] = c_tol
    else
        @error "NLP solver not setup"
    end

    x = MOI.add_variables(solver, prob.n_nlp)

    for i = 1:prob.n_nlp
        xi = MOI.SingleVariable(x[i])
        MOI.add_constraint(solver, xi, MOI.LessThan(x_u[i]))
        MOI.add_constraint(solver, xi, MOI.GreaterThan(x_l[i]))
        MOI.set(solver, MOI.VariablePrimalStart(), x[i], x0[i])
    end

    # Solve the problem
    MOI.set(solver, MOI.NLPBlock(), block_data)
    MOI.set(solver, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(solver)

    # Get the solution
    res = MOI.get(solver, MOI.VariablePrimal(), x)

    return res
end
