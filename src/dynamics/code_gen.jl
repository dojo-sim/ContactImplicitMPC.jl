"""
	generate_base_expressions(model::ContactDynamicsModel)
Generate fast base methods using Symbolics symbolic computing tools.
"""
function generate_base_expressions(model::ContactDynamicsModel;
	M_analytical = true,
	mapping = x -> I,
	nv = model.dim.q)

	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b
	np = dim(model.env)

	# Declare variables
	@variables q[1:nq]
	@variables q̇[1:nv]

	# Lagrangian
	L = lagrangian(model, q, q̇)
	L = Symbolics.simplify.(L)

	dLq = mapping(q)' * Symbolics.gradient(L, q, simplify=true) # including mapping for orientation (e.g., attitude Jacobian)
	dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
	ddL = Symbolics.hessian(L, [q; q̇], simplify=true)
	ddLq̇q = ddL[nq .+ (1:nv), 1:nq] * mapping(q)

	# Mass Matrix
	if M_analytical
		M = M_func(model, q)
		M = reshape(M, (nq, nq))
		M = Symbolics.simplify.(M)
	else
		M = ddL[nq .+ (1:nv), nq .+ (1:nv)]
	end

	# Control input Jacobian
	B = B_func(model, q)
	B = reshape(B, (nu, nv))
	B = Symbolics.simplify.(B)

	# Disturbance input Jacobian
	A = A_func(model, q)
	A = reshape(A, (nw, nv))
	A = Symbolics.simplify.(A)

	# Contact Jacobian
	J = J_func(model, q)
	J = reshape(J, (np * nc, nv))
	J = Symbolics.simplify.(J)

	# Coriolis and Centrifugal forces Jacobians
	C = ddLq̇q * q̇ - dLq
	C = Symbolics.simplify.(C)

	# Kinematics
	k = kinematics(model, q)
	k = simplify.(k)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = build_function([L], q, q̇)[1] # need to transpose to get a line vector
	expr[:M]    = build_function(M, q)[1]
	expr[:B]    = build_function(B, q)[1]
	expr[:A]    = build_function(A, q)[1]
	expr[:J]    = build_function(J, q)[1]
	expr[:C]    = build_function(C, q, q̇)[1]
	expr[:k]    = build_function(k, q)[1]

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

	# Kinematics
	k = kinematics(model, q)
	k = simplify.(k)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = :(0.0 + 0.0) # TODO: replace with base instantiation
	expr[:M]    = build_function(M, q)[1]
	expr[:B]    = build_function(B, q)[1]
	expr[:A]    = build_function(A, q)[1]
	expr[:J]    = build_function(J, q)[1]
	expr[:C]    = build_function(C, q, q̇)[1]
	expr[:k]    = build_function(k, q)[1]

	return expr
end

"""
	generate_dynamics_expressions(model::ContactDynamicsModel)
Generate fast dynamics methods using Symbolics symbolic computing tools.
"""
function generate_dynamics_expressions(model::ContactDynamicsModel; derivs = false)
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
	@variables λ1[1:ncf]
	@variables q2[1:nq]
	@variables h[1:1]

	# Expressions
	expr = Dict{Symbol, Expr}()

	# Dynamics
	d = dynamics(model, h, q0, q1, u1, w1, λ1, q2)
	d = Symbolics.simplify.(d)

	# Functions
	expr[:d] = build_function(d, h, q0, q1, u1, w1, λ1, q2)[1]

	if derivs
		dq2 = Symbolics.jacobian(d, q2, simplify = true)
		dλ1 = Symbolics.jacobian(d, λ1, simplify = true)
		dθ = Symbolics.jacobian(d, [q0; q1; u1; w1; h])#, simplify = true)

		expr[:dq2]  = build_function(dq2, h, q0, q1, u1, w1, λ1, q2)[1]
		expr[:dλ1]  = build_function(dλ1, h, q0, q1, u1, w1, λ1, q2)[1]
		expr[:dθ]   = build_function(dθ, h, q0, q1, u1, w1, λ1, q2)[1]
	end

	return expr
end

function generate_contact_expressions(model::ContactDynamicsModel;
		T = Float64, jacobians = false)

	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	nb = model.dim.b
	ncf = nc * dim(model.env)
	nz = num_var(model)
	nθ = num_data(model)

	# Declare variables
	@variables q0[1:nq]
	@variables q1[1:nq]
	@variables γ1[1:nc]
	@variables b1[1:nb]
	@variables λ1[1:ncf]
	@variables q2[1:nq]
	@variables vt[1:nb]
	@variables ψ1[1:nc]
	@variables η1[1:nb]
	@variables h[1:1]
	@variables k[1:ncf]

	@variables z[1:nz]
	@variables θ[1:nθ]
	@variables κ

	# Expressions
	expr = Dict{Symbol, Expr}()

	# Signed distance
	ϕ = ϕ_func(model, q2)
	ϕ = Symbolics.simplify.(ϕ)[1:nc]

	# Contact forces
	cf = contact_forces(model, γ1, b1, q2, k)
	cf = Symbolics.simplify.(cf)

	# Velocity stack
	vs = velocity_stack(model, q1, q2, k, h)
	vs = Symbolics.simplify.(vs)

	# Functions
	expr[:ϕ]    = build_function(ϕ, q2)[1]
	expr[:cf]   = build_function(cf, γ1, b1, q2, k)[1]
	expr[:vs]   = build_function(vs, q1, q2, k, h)[1]

	if jacobians
		dcf  = Symbolics.jacobian(cf, [q2; γ1; b1], simplify=true) # NOTE: input order change

		vsq2 = Symbolics.jacobian(vs, q2, simplify = true)
		vsq1h = Symbolics.jacobian(vs, [q1; h], simplify = true)

		# Maximum dissipation (eq.)
		md = vt + transpose(E_func(model)) * ψ1 - η1
		md = Symbolics.simplify.(md)
		mdvs = Symbolics.jacobian(md, vt, simplify = true)
		mdψη = Symbolics.jacobian(md, [ψ1; η1], simplify = true)

		# # Residual constraints
		rc = res_con(model, z, θ, κ)
		rc = Symbolics.simplify.(rc)
		rcz = Symbolics.jacobian(rc, z, simplify = true)
		rcθ = Symbolics.jacobian(rc, θ, simplify = true)

		# Functions
		expr[:dcf]  = build_function(dcf, γ1, b1, q2, k)[1]

		expr[:vsq2]  = build_function(vsq2, q1, q2, k, h)[1]
		expr[:vsq1h] = build_function(vsq1h, q1, q2, k, h)[1]

		expr[:mdvs] = build_function(mdvs, vt, ψ1, η1)[1]
		expr[:mdψη] = build_function(mdψη, vt, ψ1, η1)[1]

		expr[:rc] = build_function(rc, z, θ, κ)[2]
		expr[:rcz] = build_function(rcz, z, θ)[2]
		expr[:rcθ] = build_function(rcθ, z, θ)[2]
	end

	return expr
end

"""
	generate_residual_expressions(model::ContactDynamicsModel)
Generate fast residual methods using Symbolics symbolic computing tools.
"""
function generate_residual_expressions(model::ContactDynamicsModel;
		jacobians = :full, T = Float64)

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

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:r]  = build_function(r, z, θ, κ)[2]

	if jacobians == :full
		# contact expressions
		expr_contact = generate_contact_expressions(model,
			T = T, jacobians = false)

		rz = Symbolics.jacobian(r, z, simplify = true)
		rθ = Symbolics.jacobian(r, θ, simplify = true) # TODO: sparse version

		rz_sp = similar(rz, T)
		rθ_sp = similar(rθ, T)

		expr[:rz] = build_function(rz, z, θ)[2]
		expr[:rθ] = build_function(rθ, z, θ)[2]
	else
		expr_contact = generate_contact_expressions(model,
			T = T, jacobians = true)

		rz_sp = zeros(nz, nz)
		rθ_sp = zeros(nz, nθ)
	end

	expr = merge(expr, expr_contact)

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
	@variables   κ
    @variables  z0[1:nz]
    @variables  θ0[1:nθ]
    @variables  r0[1:nz]
    @variables rz0[1:nz,1:nz]
    @variables rθ0[1:nz,1:nθ]

	r = zeros(eltype(z), nz)
	r .= r0 + rz0 * (z-z0) + rθ0 * (θ-θ0)
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
	fct.M    = eval(expr[:M])
	fct.B    = eval(expr[:B])
	fct.A    = eval(expr[:A])
	fct.J    = eval(expr[:J])
	fct.C    = eval(expr[:C])
	fct.k    = eval(expr[:k])
	return nothing
end

"""
	instantiate_dynamics!(model,
		path::AbstractString="dynamics.jld2")
Loads the dynamics expressoins from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_dynamics!(model::ContactDynamicsModel, path::AbstractString="model/dynamics.jld2";
	derivs = false)
	expr = load_expressions(path)
	instantiate_dynamics!(model.dyn, expr, derivs = derivs)
	return nothing
end

"""
	instantiate_dynamics!(model,
		path::AbstractString="dynamics.jld2")
Evaluates the dynamics expressions to generate functions, stores them into the dedicated struture.
"""
function instantiate_dynamics!(fct::DynamicsMethods, expr::Dict{Symbol,Expr};
	derivs = false)
	fct.d = eval(expr[:d])

	if derivs
		fct.dθ = eval(expr[:dθ])
		fct.dλ1 = eval(expr[:dλ1])
		fct.dq2 = eval(expr[:dq2])
	end

	return nothing
end

"""
	instantiate_residual!(model::QuadrupedModel,
		path::AbstractString=".expr/quadruped_residual.jld2")
Loads the residual expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_residual!(model::ContactDynamicsModel, path::AbstractString="model/residual.jld2";
	jacobians = :full)
	expr = load_expressions(path)
	instantiate_contact_methods!(model.con, expr, jacobians = jacobians)
	# instantiate_residual!(model.res, expr, jacobians = jacobians)

	model.res.r!  = eval(expr[:r])

	if jacobians == :full
		model.res.rz! = eval(expr[:rz])
		model.res.rθ! = eval(expr[:rθ])
	else
		_rz(rz, z, θ) = rz_approx!(model, rz, z, θ)
		_rθ(rθ, z, θ) = rθ_approx!(model, rθ, z, θ)
		model.res.rz! = _rz
		model.res.rθ! = _rθ
	end
	return nothing
end

"""
	instantiate_residual!(model,
		path::AbstractString="model/residual.jld2")
Evaluates the residual expressions to generate functions, stores them into the model.
"""
function instantiate_residual!(fct::ResidualMethods, expr::Dict{Symbol,Expr};
	jacobians = :full)

	fct.r!  = eval(expr[:r])

	if jacobians == :full
		fct.rz! = eval(expr[:rz])
		fct.rθ! = eval(expr[:rθ])
	else
		# fct.rz! = rz_approx!
		# fct.rθ! = rθ_approx!
	end

	return nothing
end

function instantiate_contact_methods!(fct::ContactMethods, expr::Dict{Symbol,Expr};
	jacobians = :full)

	fct.ϕ = eval(expr[:ϕ])
	fct.cf = eval(expr[:cf])
	fct.vs = eval(expr[:vs])

	if jacobians == :approx
		fct.dcf = eval(expr[:dcf])
		fct.vsq2 = eval(expr[:vsq2])
		fct.vsq1h = eval(expr[:vsq1h])
		fct.mdvs = eval(expr[:mdvs])
		fct.mdψη = eval(expr[:mdψη])
		fct.rc = eval(expr[:rc])
		fct.rcz = eval(expr[:rcz])
		fct.rcθ = eval(expr[:rcθ])
	end

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
