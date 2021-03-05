"""
	generate_dynamics_expressions(model::ContactDynamicsModel)
Generate fast dynamics methods using ModelingToolkit symbolic computing tools.
"""
function generate_dynamics_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b

	# Declare variables
	@variables q0[1:nq]
	@variables q1[1:nq]
	@variables u1[1:nu]
	@variables w1[1:nw]
	@variables γ1[1:nc]
	@variables b1[1:nb]
	@variables q2[1:nq]
	@variables h[1:1]

	# Dynamics
	d = dynamics(model, h[1], q0, q1, u1, w1, γ1, b1, q2)
	d = ModelingToolkit.simplify.(d)
	dy  = ModelingToolkit.jacobian(d, [q0; q1; u1; w1; γ1; b1; q2], simplify=true)
	dq0 = ModelingToolkit.jacobian(d, q0, simplify=true)
	dq1 = ModelingToolkit.jacobian(d, q1,  simplify=true)
	du1 = ModelingToolkit.jacobian(d, u1,  simplify=true)
	dγ1 = ModelingToolkit.jacobian(d, γ1,  simplify=true)
	db1 = ModelingToolkit.jacobian(d, b1,  simplify=true)
	dq2 = ModelingToolkit.jacobian(d, q2,  simplify=true)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:d]   = build_function(d,   h, q0, q1, u1, w1, γ1, b1, q2)[1]
	expr[:dy]  = build_function(dy,  h, q0, q1, u1, w1, γ1, b1, q2)[2]
	expr[:dq0] = build_function(dq0, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	expr[:dq1] = build_function(dq1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	expr[:du1] = build_function(du1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	expr[:dγ1] = build_function(dγ1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	expr[:db1] = build_function(db1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	expr[:dq2] = build_function(dq2, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	return expr
end


"""
	generate_base_expressions(model::ContactDynamicsModel)
Generate fast base methods using ModelingToolkit symbolic computing tools.
"""
function generate_base_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b

	# Declare variables
	@variables q[1:nq]
	@variables q̇[1:nq]

	# Lagrangian
	L = lagrangian(model, q, q̇)
	L = ModelingToolkit.simplify.(L)

	dLq = ModelingToolkit.gradient(L, q, simplify=true)
	dLq̇ = ModelingToolkit.gradient(L, q̇, simplify=true)
	ddL = ModelingToolkit.sparsehessian(L, [q; q̇], simplify=true)
	ddL = SparseMatrixCSC{Expression,Int64}(ddL)
	ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

	# Mass Matrix
	M = M_func(model, q)
	M = ModelingToolkit.simplify.(M)

	# Control input Jacobian
	B = B_func(model, q)
	B = reshape(B, (nb, nq))
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
function generate_residual_expressions(model::ContactDynamicsModel; T = Float64)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nz = num_var(model)
	nθ = num_data(model)

	# Declare variables
	@variables z[1:nz]
	@variables θ[1:nθ]
	@variables κ[1:1]

	# Residual
	r = residual(model, z, θ, κ[1])
	r = ModelingToolkit.simplify.(r)
	rz = ModelingToolkit.sparsejacobian(r, z, simplify = true)
	rθ = ModelingToolkit.jacobian(r, θ, simplify = true) # TODO: sparse version

	rz_sp = similar(rz, T)
	rθ_sp = similar(rθ, T)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:r]  = build_function(r, z, θ, κ)[2]
	expr[:rz] = build_function(rz, z, θ, κ)[2]
	expr[:rθ] = build_function(rθ, z, θ, κ)[2]
	return expr, rz_sp, rθ_sp
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

function instantiate_base(expr::Dict{Symbol,Expr})
	return eval(expr[:L]), eval(expr[:M]), eval(expr[:B]), eval(expr[:N]), eval(expr[:P]), eval(expr[:C])
end

function instantiate_dynamics(expr::Dict{Symbol,Expr})
	eval(expr[:d]), eval(expr[:dy]), eval(expr[:dq0]), eval(expr[:dq1]), eval(expr[:du1]), eval(expr[:dγ1]), eval(expr[:db1]), eval(expr[:dq2])
end

function instantiate_residual(expr::Dict{Symbol,Expr})
	eval(expr[:r]), eval(expr[:rz]), eval(expr[:rθ])
end

function generate_fast_expressions(model::ContactDynamicsModel, name::String;
		generate = true, verbose = true)
	model = deepcopy(model)
	verbose && println()
	if generate
		verbose && println("generating base methods...")
		expr_base = generate_base_expressions(model)
		save_expressions(expr_base, joinpath(@__DIR__, name, "base.jld2"), overwrite=true)

		verbose && println("generating dynamics methods...")
		expr_dyn = generate_dynamics_expressions(model)
		save_expressions(expr_dyn, joinpath(@__DIR__, name, "dynamics.jld2"), overwrite=true)

		verbose && println("generating residual methods...")
		expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
		save_expressions(expr_res, joinpath(@__DIR__, name, "residual.jld2"), overwrite=true)
		@save joinpath(@__DIR__, name, "sparse_jacobians.jld2") rz_sp rθ_sp
	else
		verbose && println("loading base methods...")
		expr_base = load_expressions(joinpath(@__DIR__, name, "base.jld2"))

		verbose && println("loading dynamics methods...")
		expr_dyn = load_expressions(joinpath(@__DIR__, name, "dynamics.jld2"))

		verbose && println("loading residual methods...")
		expr_res = load_expressions(joinpath(@__DIR__, name, "residual.jld2"))
		@load joinpath(@__DIR__, name, "sparse_jacobians.jld2") rz_sp rθ_sp
	end

	verbose && println("generated methods: success")

	L_fast, M_fast, B_fast, N_fast, P_fast, C_fast = instantiate_base(expr_base)
	d, dy!, dq0!, dq1!, du1!, dγ1!, db1!, dq2! = instantiate_dynamics(expr_dyn)
	r!, rz!, rθ! = instantiate_residual(expr_res)

	verbose && println("instantiated methods: success")

	return L_fast, M_fast, B_fast, N_fast, P_fast, C_fast, d, dy!, dq0!, dq1!, du1!, dγ1!, db1!, dq2!, r!, rz!, rθ!, rz_sp, rθ_sp
end
