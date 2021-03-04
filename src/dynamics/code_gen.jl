
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
	@variables q0[1:nq]
	@variables q1[1:nq]
	@variables u1[1:nu]
	@variables γ1[1:nγ]
	@variables b1[1:nb]
	@variables q2[1:nq]
	@variables dt[1:1]

	# Dynamics
	d = dynamics(model,dt[1],q0,q1,u1,γ1,b1,q2)
	d = ModelingToolkit.simplify.(d)
	dy  = ModelingToolkit.jacobian(d, [q0; q1; u1; γ1; b1; q2], simplify=true)
	dq0 = ModelingToolkit.jacobian(d, q0, simplify=true)
	dq1 = ModelingToolkit.jacobian(d, q1,   simplify=true)
	du1 = ModelingToolkit.jacobian(d, u1,   simplify=true)
	dγ1 = ModelingToolkit.jacobian(d, γ1,   simplify=true)
	db1 = ModelingToolkit.jacobian(d, b1,   simplify=true)
	dq2 = ModelingToolkit.jacobian(d, q2,  simplify=true)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:d]   = build_function(d,   dt, q0, q1, u1, γ1, b1, q2)[1]
	expr[:dy]  = build_function(dy,  dt, q0, q1, u1, γ1, b1, q2)[2]
	expr[:dq0] = build_function(dq0, dt, q0, q1, u1, γ1, b1, q2)[2]
	expr[:dq1] = build_function(dq1, dt, q0, q1, u1, γ1, b1, q2)[2]
	expr[:du1] = build_function(du1, dt, q0, q1, u1, γ1, b1, q2)[2]
	expr[:dγ1] = build_function(dγ1, dt, q0, q1, u1, γ1, b1, q2)[2]
	expr[:db1] = build_function(db1, dt, q0, q1, u1, γ1, b1, q2)[2]
	expr[:dq2] = build_function(dq2, dt, q0, q1, u1, γ1, b1, q2)[2]
	return expr
end


"""
	generate_base_expressions(model::ContactDynamicsModel)
Generate fast base methods using ModelingToolkit symbolic computing tools.
"""
function generate_base_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b
	# Declare variables
	@variables q[1:nq]
	@variables q̇[1:nq]
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
	B = reshape(B, (nb,nq))
	B = ModelingToolkit.simplify.(B)

	# Impact force Jacobian
	N = N_func(model, q)
	N = ModelingToolkit.simplify.(N)

	# Friction Force Jacobian
	P = P_func(model, q)
	P = ModelingToolkit.simplify.(P)

	# Coriolis and Centrifugal forces Jacobians
	C = ddLq̇q * q̇ - dLq
	C = ModelingToolkit.simplify.(C)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = build_function(transpose([L]), q, q̇)[1] # need to transpose to get a line vector
	expr[:M]    = build_function(M, q)[1]
	expr[:B]    = build_function(B, q)[1]
	expr[:N]    = build_function(N, q)[1]
	expr[:P]    = build_function(P, q)[1]
	expr[:C]    = build_function(C, q, q̇)[1]
	return expr
end


"""
	generate_residual_expressions(model::ContactDynamicsModel)
Generate fast residual methods using ModelingToolkit symbolic computing tools.
"""
function generate_residual_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nγ = model.dim.γ
	nb = model.dim.b
	nz = model.dim.z
	nθ = model.dim.θ

	# Declare variables
	@variables dt[1:1]
	@variables z[1:nz]
	@variables θ[1:nθ]
	@variables ρ[1:1]
	# Residual
	r = residual(model, dt[1], z, θ, ρ[1])
	r = ModelingToolkit.simplify.(r)
	rz = ModelingToolkit.jacobian(r, z, simplify=true)
	rθ = ModelingToolkit.jacobian(r, θ, simplify=true)


	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:r]  = build_function(r,  dt, z, θ, ρ)[2]
	expr[:rz] = build_function(rz, dt, z, θ, ρ)[2]
	expr[:rθ] = build_function(rθ, dt, z, θ, ρ)[2]
	return expr
end


"""
	save_expressions(expr::Dict{Symbol,Expr},
		path::AbstractString="dynamics_expressions.jld2"; overwrite::Bool=false)
Save the fast expressions obtained with ModelingToolkit in a `jld2` file.
"""
function save_expressions(expr::Dict{Symbol,Expr},
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
	load_expressions(path::AbstractString="dynamics_expressions.jld2")
Load the fast expressions obtained with ModelingToolkit from a `jld2` file.
"""
function load_expressions(path::AbstractString="dynamics_expressions.jld2")
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
function instantiate_dynamics!(model::QuadrupedModel, path::AbstractString=".expr/quadruped_dynamics.jld2")
	expr = load_expressions(path)
	instantiate_dynamics!(model.dyn, expr)
	return nothing
end

"""
	instantiate_dynamics!(model::QuadrupedModel,
		path::AbstractString="quadruped_expr.jld2")
Evaluates the dynamics expressions to generate functions, stores them into the dedicated struture.
"""
function instantiate_dynamics!(fct::DynamicsMethods13, expr::Dict{Symbol,Expr})
	fct.d   = eval(expr[:d])
	fct.dy  = eval(expr[:dy])
	fct.dq0 = eval(expr[:dq0])
	fct.dq1 = eval(expr[:dq1])
	fct.du1 = eval(expr[:du1])
	fct.dγ1 = eval(expr[:dγ1])
	fct.db1 = eval(expr[:db1])
	fct.dq2 = eval(expr[:dq2])
	return nothing
end

"""
	instantiate_base!(model::QuadrupedModel,
		path::AbstractString=".expr/quadruped_base.jld2")
Loads the base expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_base!(model::QuadrupedModel, path::AbstractString=".expr/quadruped_base.jld2")
	expr = load_expressions(path)
	instantiate_base!(model.bas, expr)
	return nothing
end

"""
	instantiate_base!(model::QuadrupedModel,
		path::AbstractString=".expr/quadruped_base.jld2")
Evaluates the base expressions to generate functions, stores them into the model.
"""
function instantiate_base!(fct::BaseMethods12, expr::Dict{Symbol,Expr})
	fct.L    = eval(expr[:L])
	fct.M    = eval(expr[:M])
	fct.B    = eval(expr[:B])
	fct.N    = eval(expr[:N])
	fct.P    = eval(expr[:P])
	fct.C    = eval(expr[:C])
	return nothing
end

"""
	instantiate_residual!(model::QuadrupedModel,
		path::AbstractString=".expr/quadruped_residual.jld2")
Loads the residual expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_residual!(model::QuadrupedModel, path::AbstractString=".expr/quadruped_residual.jld2")
	expr = load_expressions(path)
	instantiate_residual!(model.res, expr)
	return nothing
end

"""
	instantiate_residual!(model::QuadrupedModel,
		path::AbstractString=".expr/quadruped_residual.jld2")
Evaluates the residual expressions to generate functions, stores them into the model.
"""
function instantiate_residual!(fct::ResidualMethods14, expr::Dict{Symbol,Expr})
	fct.r   = eval(expr[:r])
	fct.rz = eval(expr[:rz])
	fct.rθ = eval(expr[:rθ])
	return nothing
end
