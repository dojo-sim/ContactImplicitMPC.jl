
"""
	generate_dynamics_expressions(model::ContactDynamicsModel)
Generate fast dynamics methods using ModelingToolkit symbolic computing tools.
"""
function generate_dynamics_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b
	# Declare variables
	@variables q_1[1:nq]
	@variables q[1:nq]
	@variables u[1:nu]
	@variables γ[1:nγ]
	@variables b[1:nb]
	@variables q1[1:nq]
	@variables q̇[1:nq]
	@variables dt[1:1]

	# Symbolic computing
	# Lagrangian
	L = lagrangian(model, q, q̇)
	L = ModelingToolkit.simplify.(L)

	dLq = ModelingToolkit.gradient(L, q, simplify=true)
	dLq̇ = ModelingToolkit.gradient(L, q̇, simplify=true)
	ddL = ModelingToolkit.sparsehessian(L, [q;q̇], simplify=true)
	ddL = SparseMatrixCSC{Expression,Int64}(ddL)
	ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

	# Mass Matrix
	M = M_func(model, q)
	M = ModelingToolkit.simplify.(M)

	# Control input Jacobian
	B = B_func(model, q)
	B = ModelingToolkit.simplify.(B)

	# Impact force Jacobian
	N = N_func(model, q)
	N = ModelingToolkit.simplify.(N)

	# Friction Force Jacobian
	P = _P_func(model, q)
	P = ModelingToolkit.simplify.(P)

	# Coriolis and Centrifugal forces Jacobians
	C = ddLq̇q * q̇ - dLq
	C = ModelingToolkit.simplify.(C)

	# Dynamics
	d = dynamics(model,dt[1],q_1,q,u,γ,b,q1)
	d = ModelingToolkit.simplify.(d)
	dz   = ModelingToolkit.jacobian(d, [q_1; q; u; γ; b; q1], simplify=true)
	dq_1 = ModelingToolkit.jacobian(d, q_1, simplify=true)
	dq   = ModelingToolkit.jacobian(d, q,   simplify=true)
	du   = ModelingToolkit.jacobian(d, q,   simplify=true)
	dγ   = ModelingToolkit.jacobian(d, γ,   simplify=true)
	db   = ModelingToolkit.jacobian(d, b,   simplify=true)
	dq1  = ModelingToolkit.jacobian(d, q1,  simplify=true)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = build_function(transpose([L]), q, q̇)[1] # need to transpose to get a line vector
	expr[:M]    = build_function(M, q)[1]
	expr[:B]    = build_function(B, q)[1]
	expr[:N]    = build_function(N, q)[1]
	expr[:P]    = build_function(P, q)[1]
	expr[:C]    = build_function(C, q, q̇)[1]
	expr[:d]    = build_function(d,    dt, q_1, q, u, γ, b, q1)[1]
	expr[:dz]   = build_function(dz,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:dq_1] = build_function(dq_1, dt, q_1, q, u, γ, b, q1)[2]
	expr[:dq]   = build_function(dq,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:du]   = build_function(du,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:dγ]   = build_function(dγ,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:db]   = build_function(db,   dt, q_1, q, u, γ, b, q1)[2]
	expr[:dq1]  = build_function(dq1,  dt, q_1, q, u, γ, b, q1)[2]
	return expr
end


"""
	save_dynamics_expressions(expr::Dict{Symbol,Expr},
		path::AbstractString="dynamics_expressions.jld2"; overwrite::Bool=false)
Save the fast dynamics expressions obtained with ModelingToolkit in a `jld2` file.
"""
function save_dynamics_expressions(expr::Dict{Symbol,Expr},
	path::AbstractString="dynamics_expressions.jld2"; overwrite::Bool=false)
	path = abspath(path)
    if isfile(path) && !overwrite
        error("The output path $path already exists. To overwrite that file, you can pass `overwrite=true` to this function")
    end
	@save path expr
	@info("Saved output as $path")
	return nothing
end

"""
	load_dynamics_expressions(path::AbstractString="dynamics_expressions.jld2")
Load the fast dynamics expressions obtained with ModelingToolkit from a `jld2` file.
"""
function load_dynamics_expressions(path::AbstractString="dynamics_expressions.jld2")
	path = abspath(path)
	if !isfile(path)
		error("Could not find the input file $path")
	end
	@load path expr
	return expr
end

"""
	instantiate_dynamics!(model::QuadrupedModel,
		path::AbstractString="quadruped_expr.jld2")
Loads the dynamics expressoins from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_dynamics!(model::QuadrupedModel, path::AbstractString="quadruped_expr.jld2")
	expr = load_dynamics_expressions(path)
	instantiate_dynamics!(model.fct, expr)
	return nothing
end

"""
	instantiate_dynamics!(model::QuadrupedModel,
		path::AbstractString="quadruped_expr.jld2")
Evaluates the dynamics expressions to generate functions, stores them into the dedicated struture.
"""
function instantiate_dynamics!(fct::DynamicsMethods11, expr::Dict{Symbol,Expr})
	fct.L    = eval(expr[:L])
	fct.M    = eval(expr[:M])
	fct.B    = eval(expr[:B])
	fct.N    = eval(expr[:N])
	fct.P    = eval(expr[:P])
	fct.C    = eval(expr[:C])
	fct.d    = eval(expr[:d])
	fct.dz   = eval(expr[:dz])
	fct.dq_1 = eval(expr[:dq_1])
	fct.dq   = eval(expr[:dq])
	fct.du   = eval(expr[:du])
	fct.dγ   = eval(expr[:dγ])
	fct.db   = eval(expr[:db])
	fct.dq1  = eval(expr[:dq1])
	return nothing
end
