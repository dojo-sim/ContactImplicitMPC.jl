"""
	generate_base_expressions(model::ContactDynamicsModel)
Generate fast base methods using Symbolics symbolic computing tools.
"""
function generate_base_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b
	np = dim(model.env)

	# Declare variables
	@variables q[1:nq]
	@variables q̇[1:nq]

	# Lagrangian
	L = lagrangian(model, q, q̇)
	L = Symbolics.simplify.(L)

	dLq = Symbolics.gradient(L, q, simplify=true)
	dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
	ddL = Symbolics.sparsehessian(L, [q; q̇], simplify=true)
	ddLq̇q = ddL[nq .+ (1:nq), 1:nq]

	# Signed distance
	ϕ = ϕ_func(model, q)
	ϕ = Symbolics.simplify.(ϕ)[1:nc]

	# Mass Matrix
	M = M_func(model, q)
	M = reshape(M, nq, nq)
	M = Symbolics.simplify.(M)

	# Control input Jacobian
	B = B_func(model, q)
	B = reshape(B, (nu, nq))
	B = Symbolics.simplify.(B)

	# Control input Jacobian
	A = A_func(model, q)
	A = reshape(A, (nw, nq))
	A = Symbolics.simplify.(A)

	# Contact Force Jacobian
	J = J_func(model, q)
	J = reshape(J, np * nc, nq)
	J = Symbolics.simplify.(J)

	# Coriolis and Centrifugal forces Jacobians
	C = ddLq̇q * q̇ - dLq
	C = Symbolics.simplify.(C)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = build_function([L], q, q̇)[1] # need to transpose to get a line vector
	expr[:ϕ]    = build_function(ϕ, q)[1]
	expr[:M]    = build_function(M, q)[1]
	expr[:B]    = build_function(B, q)[1]
	expr[:A]    = build_function(A, q)[1]
	expr[:J]    = build_function(J, q)[1]
	expr[:C]    = build_function(C, q, q̇)[1]
	return expr
end

function generate_base_expressions_analytical(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b
	np = dim(model.env)

	# Declare variables
	@variables q[1:nq]
	@variables q̇[1:nq]

	# Signed distance
	ϕ = ϕ_func(model, q)
	ϕ = Symbolics.simplify.(ϕ)[1:nc]

	# Mass Matrix
	M = M_func(model, q)
	M = reshape(M, nq, nq)
	M = Symbolics.simplify.(M)

	# Control input Jacobian
	B = B_func(model, q)
	B = reshape(B, (nu, nq))
	B = Symbolics.simplify.(B)

	# Control input Jacobian
	A = A_func(model, q)
	A = reshape(A, (nw, nq))
	A = Symbolics.simplify.(A)

	# Contact Force Jacobian
	J = J_func(model, q)
	J = reshape(J, np * nc, nq)
	J = Symbolics.simplify.(J)

	# Coriolis and Centrifugal forces Jacobians
	C = C_func(model, q, q̇)
	C = Symbolics.simplify.(C)[1:nq]

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = :(0.0 + 0.0) # TODO: replace with base instantiation
	expr[:ϕ]    = build_function(ϕ, q)[1]
	expr[:M]    = build_function(M, q)[1]
	expr[:B]    = build_function(B, q)[1]
	expr[:A]    = build_function(A, q)[1]
	expr[:J]    = build_function(J, q)[1]
	expr[:C]    = build_function(C, q, q̇)[1]
	return expr
end

"""
	generate_dynamics_expressions(model::ContactDynamicsModel)
Generate fast dynamics methods using Symbolics symbolic computing tools.
"""
function generate_dynamics_expressions(model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b
	ncf = nc * dim(model.env)

	# Declare variables
	@variables q0[1:nq]
	@variables q1[1:nq]
	@variables u1[1:nu]
	@variables w1[1:nw]
	@variables γ1[1:nc]
	@variables λ1[1:ncf]
	@variables b1[1:nb]
	@variables q2[1:nq]
	@variables h

	# Dynamics
	# d = dynamics(model, h, q0, q1, u1, w1, γ1, b1, q2)
	d = dynamics(model, h, q0, q1, u1, w1, λ1, q2)
	d = Symbolics.simplify.(d)
	# dy  = Symbolics.jacobian(d, [q0; q1; u1; w1; γ1; b1; q2], simplify=true)
	# dq0 = Symbolics.jacobian(d, q0, simplify=true)
	# dq1 = Symbolics.jacobian(d, q1,  simplify=true)
	# du1 = Symbolics.jacobian(d, u1,  simplify=true)
	# dw1 = Symbolics.jacobian(d, w1,  simplify=true)
	# dγ1 = Symbolics.jacobian(d, γ1,  simplify=true)
	# db1 = Symbolics.jacobian(d, b1,  simplify=true)
	# dq2 = Symbolics.jacobian(d, q2,  simplify=true)

	# Build function
	expr = Dict{Symbol, Expr}()
	# expr[:d]   = build_function(d,   h, q0, q1, u1, w1, γ1, b1, q2)[1]
	expr[:d]   = build_function(d,   h, q0, q1, u1, w1, λ1, q2)[1]
	# expr[:dy]  = build_function(dy,  h, q0, q1, u1, w1, γ1, b1, q2)[2]
	# expr[:dq0] = build_function(dq0, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	# expr[:dq1] = build_function(dq1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	# expr[:du1] = build_function(du1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	# expr[:dw1] = build_function(dw1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	# expr[:dγ1] = build_function(dγ1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	# expr[:db1] = build_function(db1, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	# expr[:dq2] = build_function(dq2, h, q0, q1, u1, w1, γ1, b1, q2)[2]
	return expr
end

"""
	generate_residual_expressions(model::ContactDynamicsModel)
Generate fast residual methods using Symbolics symbolic computing tools.
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
	@variables κ

	# Residual
	r = residual(model, z, θ, κ)
	r = Symbolics.simplify.(r)
	rz = Symbolics.jacobian(r, z, simplify = true)
	rθ = Symbolics.jacobian(r, θ, simplify = true) # TODO: sparse version

	rz_sp = similar(rz, T)
	rθ_sp = similar(rθ, T)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:r]  = build_function(r, z, θ, κ)[2]
	expr[:rz] = build_function(rz, z, θ)[2]
	expr[:rθ] = build_function(rθ, z, θ)[2]
	return expr, rz_sp, rθ_sp
end

"""
	generate_linearized_expressions(model::ContactDynamicsModel)
Generate fast linearizedimate methods using Symbolics symbolic computing tools.
"""
function generate_linearized_expressions(model::ContactDynamicsModel; T = Float64)
    bil_terms, bil_vars = get_bilinear_indices(model)
    nz = num_var(model)
    nθ = num_data(model)

	# r_linearized rz_linearized, rθ_linearized
    @variables   z[1:nz]
	@variables   θ[1:nθ]
	@variables   κ[1:1]
    @variables  z0[1:nz]
    @variables  θ0[1:nθ]
    @variables  r0[1:nz]
    @variables rz0[1:nz,1:nz]
    @variables rθ0[1:nz,1:nθ]

	r = zeros(eltype(z), nz)
	r .= r0 + rz0 * (z-z0) + rθ0 * (θ-θ0)
    # r .= rz0 * (z-z0) + rθ0 * (θ-θ0) # wrong
	for i = 1:length(bil_terms)
		t = bil_terms[i]
		v1 = bil_vars[i][1]
		v2 = bil_vars[i][2]
		for j = 1:length(t)
			r[t[j]] = z[v1[j]]*z[v2[j]] - κ[1]
		end
	end
	r = simplify.(r)

	# r_linearized rz_linearized, rθ_linearized
    @variables   z[1:nz]
    @variables rz0[1:nz,1:nz]

	rz = zeros(eltype(z), nz, nz)
	rz .= rz0
	for i = 1:length(bil_terms)
		t = bil_terms[i]
		v1 = bil_vars[i][1]
		v2 = bil_vars[i][2]
		for j = 1:length(t)
			rz[t[j], v1[j]] = z[v2[j]]
			rz[t[j], v2[j]] = z[v1[j]]
		end
	end
	rz = simplify.(rz)

	# r_linearized rz_linearized, rθ_linearized
    @variables rθ0[1:nz,1:nθ]
	rθ = zeros(eltype(rθ0), nz, nθ)
	rθ .= rθ0
    rθ = simplify.(rθ)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:r]  = build_function(r,  z, θ, κ, z0, θ0, r0, rz0, rθ0)[2]
	expr[:rz] = build_function(rz, z, rz0)[2]
	expr[:rθ] = build_function(rθ, rθ0)[2]
	return expr
end

"""
	save_expressions(expr::Dict{Symbol,Expr},
		path::AbstractString="dynamics_expressions.jld2"; overwrite::Bool=false)
Save the fast expressions obtained with Symbolics in a `jld2` file.
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
Load the fast expressions obtained with Symbolics from a `jld2` file.
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
	instantiate_base!(model,
		path::AbstractString="model/base.jld2")
Loads the base expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_base!(model::ContactDynamicsModel, path::AbstractString="model/base.jld2")
	expr = load_expressions(path)
	instantiate_base!(model.base, expr)
	return nothing
end

"""
	instantiate_base!(model,
		path::AbstractString="model/base.jld2")
Evaluates the base expressions to generate functions, stores them into the model.
"""
function instantiate_base!(fct::BaseMethods, expr::Dict{Symbol,Expr})
	fct.L    = eval(expr[:L])
	fct.ϕ    = eval(expr[:ϕ])
	fct.M    = eval(expr[:M])
	fct.B    = eval(expr[:B])
	fct.A    = eval(expr[:A])
	fct.J    = eval(expr[:J])
	fct.C    = eval(expr[:C])
	return nothing
end

"""
	instantiate_dynamics!(model,
		path::AbstractString="dynamics.jld2")
Loads the dynamics expressoins from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_dynamics!(model::ContactDynamicsModel, path::AbstractString="model/dynamics.jld2")
	expr = load_expressions(path)
	instantiate_dynamics!(model.dyn, expr)
	return nothing
end

"""
	instantiate_dynamics!(model,
		path::AbstractString="dynamics.jld2")
Evaluates the dynamics expressions to generate functions, stores them into the dedicated struture.
"""
function instantiate_dynamics!(fct::DynamicsMethods, expr::Dict{Symbol,Expr})
	fct.d   = eval(expr[:d])
	# fct.dy  = eval(expr[:dy])
	# fct.dq0 = eval(expr[:dq0])
	# fct.dq1 = eval(expr[:dq1])
	# fct.du1 = eval(expr[:du1])
	# fct.dw1 = eval(expr[:dw1])
	# fct.dγ1 = eval(expr[:dγ1])
	# fct.db1 = eval(expr[:db1])
	# fct.dq2 = eval(expr[:dq2])
	return nothing
end

"""
	instantiate_residual!(model::QuadrupedModel,
		path::AbstractString=".expr/quadruped_residual.jld2")
Loads the residual expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_residual!(model::ContactDynamicsModel, path::AbstractString="model/residual.jld2")
	expr = load_expressions(path)
	instantiate_residual!(model.res, expr)
	return nothing
end

"""
	instantiate_residual!(model,
		path::AbstractString="model/residual.jld2")
Evaluates the residual expressions to generate functions, stores them into the model.
"""
function instantiate_residual!(fct::ResidualMethods, expr::Dict{Symbol,Expr})
	fct.r!  = eval(expr[:r])
	fct.rz! = eval(expr[:rz])
	fct.rθ! = eval(expr[:rθ])
	return nothing
end

"""
	instantiate_linearized!(model::QuadrupedModel,
		path::AbstractString=".expr/quadruped_linearized.jld2")
Loads the linearizedimate expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_linearized!(model::ContactDynamicsModel, path::AbstractString="model/linearized.jld2")
	expr = load_expressions(path)
	instantiate_residual!(model.linearized, expr)
	return nothing
end

function fast_expressions!(model::ContactDynamicsModel, dir::String;
	generate = true, verbose = true)
	verbose && println()

	model = deepcopy(model)

	path_base = joinpath(dir, "base.jld2")
	path_dyn = joinpath(dir, "dynamics.jld2")
	path_res = joinpath(dir, "residual.jld2")
	path_jac = joinpath(dir, "sparse_jacobians.jld2")

	if generate
		expr_base = generate_base_expressions(model)
		save_expressions(expr_base, path_base, overwrite=true)
		instantiate_base!(model, path_base)

		expr_dyn = generate_dynamics_expressions(model)
		save_expressions(expr_dyn, path_dyn, overwrite=true)
		instantiate_dynamics!(model, path_dyn)

		expr_res, rz_sp, rθ_sp = generate_residual_expressions(model)
		save_expressions(expr_res, path_res, overwrite=true)
		@save path_jac rz_sp rθ_sp
		instantiate_residual!(model, path_res)
		verbose && println("generated methods: success")
	else
		instantiate_base!(model, path_base)
		instantiate_dynamics!(model, path_dyn)
		instantiate_residual!(model, path_res)
		@load joinpath(dir, "sparse_jacobians.jld2") rz_sp rθ_sp
	end
	verbose && println("instantiating methods: success")
end
