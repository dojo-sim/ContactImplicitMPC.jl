
function generate_contact_expressions(model::Model, env::Environment;
		T = Float64, jacobians = false, nv = model.nq)

	# model dimensions
	nq = model.nq
	nu = model.nu
	nw = model.nw
	nc = model.nc

	# friction dimensions
	nb = nc * friction_dim(env)
	ncf = nc * dim(env)

	# simulator dimensions
	nz = num_var(model, env)
	nθ = num_data(model)

	# Declare variables
	@variables q0[1:nq]
	@variables q1[1:nq]
	@variables γ1[1:nc]
	@variables b1[1:nb]
	@variables λ1[1:ncf]
	@variables u1[1:nu]
	@variables w1[1:nw]
	@variables q2[1:nq]
	@variables vt[1:nb]
	@variables ψ1[1:nc]
	@variables η1[1:nb]
	@variables h[1:1]
	@variables k[1:ncf]

	@variables z[1:nz]
	@variables θ[1:nθ]
	@variables κ[1:1]

	# Expressions
	expr = Dict{Symbol, Expr}()

	# Contact Jacobian
	J = J_func(model, env, q2)
	J = reshape(J, (ncf, nv))
	J = Symbolics.simplify.(J)

	# Signed distance
	ϕ = ϕ_func(model, env, q2)
	ϕ = Symbolics.simplify.(ϕ)[1:nc]

	# Contact forces
	cf = contact_forces(model, env, γ1, b1, q2, k)
	cf = Symbolics.simplify.(cf)

	# Velocity stack
	vs = velocity_stack(model, env, q1, q2, k, h)
	vs = Symbolics.simplify.(vs)

	# Dynamics
	d = dynamics(model, env, h, q0, q1, u1, w1, λ1, q2)
	d = Symbolics.simplify.(d)

	# Functions
	expr[:J]    = build_function(J, q2, checkbounds=true)[1]
	expr[:ϕ]    = build_function(ϕ, q2, checkbounds=true)[1]
	expr[:cf]   = build_function(cf, γ1, b1, q2, k, checkbounds=true)[1]
	expr[:vs]   = build_function(vs, q1, q2, k, h, checkbounds=true)[1]
	expr[:d]    = build_function(d, h, q0, q1, u1, w1, λ1, q2, checkbounds=true)[1]

	if jacobians
		dJ = Symbolics.jacobian(J'*λ1, q2, simplify = true) # dq2 = ∂q2  + dJ(λ1, q2)
		dλ1 = J # this is an analytical result, might be wrong but should be okay

		dcf  = Symbolics.jacobian(cf, [q2; γ1; b1], simplify=true) # NOTE: input order change

		vsq2 = Symbolics.jacobian(-vs, q2, simplify = true)
		vsq1h = Symbolics.jacobian(-vs, [q1; h], simplify = true)

		# Maximum dissipation (eq.)
		ψ_stack = transpose(E_func(model, env)) * ψ1
		md = η1 - vt - ψ_stack
		md = Symbolics.simplify.(md)
		mdvs = Symbolics.jacobian(-md, vt, simplify = true)
		mdψη = Symbolics.jacobian(md, [ψ1; η1], simplify = true)

		# # Residual constraints
		rc = res_con(model, env, z, θ, κ)
		rc = Symbolics.simplify.(rc)
		rcz = Symbolics.jacobian(rc, z, simplify = true)
		rcθ = Symbolics.jacobian(rc, θ, simplify = true)

		# Functions
		expr[:dJ]   = build_function(dJ, λ1, q2, checkbounds=true)[1]
		expr[:dλ1]  = build_function(dλ1, q2, checkbounds=true)[1]

		expr[:dcf]  = build_function(dcf, γ1, b1, q2, k, checkbounds=true)[1]

		expr[:vsq2]  = build_function(vsq2, q1, q2, k, h, checkbounds=true)[1]
		expr[:vsq1h] = build_function(vsq1h, q1, q2, k, h, checkbounds=true)[1]

		expr[:mdvs] = build_function(mdvs, vt, ψ1, η1, checkbounds=true)[1]
		expr[:mdψη] = build_function(mdψη, vt, ψ1, η1, checkbounds=true)[1]

		expr[:rc]  = build_function(rc, z, θ, κ, checkbounds=true)[2]
		expr[:rcz] = build_function(rcz, z, θ, checkbounds=true)[2]
		expr[:rcθ] = build_function(rcθ, z, θ, checkbounds=true)[2]
	end
	return expr
end

"""
	generate_residual_expressions(model::Model)
Generate fast residual methods using Symbolics symbolic computing tools.
"""
function generate_residual_expressions(model::Model, env::Environment;
		mapping = (a,b,c) -> Diagonal(ones(num_var(model, env))), jacobians = :full,
		T = Float64, nv = model.nq)

	nq = model.nq
	nu = model.nu
	nc = model.nc
	nb = nc * friction_dim(env)
	nz = num_var(model, env)
	nθ = num_data(model)

	# Declare variables
	@variables z[1:nz]
	@variables θ[1:nθ]
	@variables κ[1:1]

	# Residual
	r = residual(model, env, z, θ, κ)
	r = Symbolics.simplify.(r)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:r]  = build_function(r, z, θ, κ, checkbounds=true)[2]

	if jacobians == :full
		# contact expressions
		expr_contact = generate_contact_expressions(model, env,
			T = T, jacobians = false, nv = nv)

		m = mapping(model, env, z)
		m = simplify.(m)
		rz = Symbolics.jacobian(r, z, simplify = true)
		rz = rz * m
		rz = simplify.(rz)


		rθ = Symbolics.jacobian(r, θ, simplify = true) # TODO: sparse version

		rz_sp = similar(rz, T)
		rθ_sp = similar(rθ, T)

		expr[:rz] = build_function(rz, z, θ, checkbounds=true)[2]
		expr[:rθ] = build_function(rθ, z, θ, checkbounds=true)[2]
	else
		expr_contact = generate_contact_expressions(model, env,
			T = T, jacobians = true, nv = nv)

		rz_sp = zeros(nz, nz)
		rθ_sp = zeros(nz, nθ)
	end

	expr = merge(expr, expr_contact)

	return expr, rz_sp, rθ_sp
end

function instantiate_residual!(s::Simulation, path_res, path_jac;
	jacobians = :full)

	# load expressions
	expr = load_expressions(path_res)
	instantiate_contact_methods!(s.con, expr, jacobians = jacobians)

	# residual
	s.res.r!  = eval(expr[:r])

	# Jacobians
	if jacobians == :full
		s.res.rz! = eval(expr[:rz])
		s.res.rθ! = eval(expr[:rθ])
	else
		_rz(rz, z, θ) = rz_approx!(s, rz, z, θ)
		_rθ(rθ, z, θ) = rθ_approx!(s, rθ, z, θ)
		s.res.rz! = _rz
		s.res.rθ! = _rθ
	end

	# pre-allocated memory
	jac = load(path_jac)
	s.rz = jac["rz_sp"]
	s.rθ = jac["rθ_sp"]

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

	fct.J = eval(expr[:J])
	fct.ϕ = eval(expr[:ϕ])
	fct.cf = eval(expr[:cf])
	fct.vs = eval(expr[:vs])
	fct.d = eval(expr[:d])

	if jacobians == :approx
		fct.dJ = eval(expr[:dJ])
		fct.dλ1 = eval(expr[:dλ1])
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
