"""
	generate_base_expressions(model::ContactModel)
Generate fast base methods using Symbolics symbolic computing tools.
"""
function generate_base_expressions(model::ContactModel;
	M_analytical = true,
	mapping = (a, x) -> I,
	nv = model.dim.q)

	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c

	# Declare variables
	@variables q[1:nq]
	@variables q̇[1:nv]

	# Lagrangian
	L = lagrangian(model, q, q̇)
	L = Symbolics.simplify.(L)

	# Mass Matrix
	if M_analytical
		M = M_func(model, q)
		M = reshape(M, (nv, nv))
		M = Symbolics.simplify.(M)

		C = C_func(model, q, q̇)
		C = simplify.(C)
	else
		dLq = mapping(model, q)' * Symbolics.gradient(L, q, simplify=true) # including mapping for orientation (e.g., attitude Jacobian)
		dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
		ddL = Symbolics.hessian(L, [q; q̇], simplify=true)
		ddLq̇q = ddL[nq .+ (1:nv), 1:nq] * mapping(model, q)

		M = ddL[nq .+ (1:nv), nq .+ (1:nv)]

		# Coriolis and Centrifugal forces Jacobians
		C = ddLq̇q * q̇ - dLq
		C = Symbolics.simplify.(C)
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
	J = reshape(J, size(J_func(model, zeros(nq))))
	J = Symbolics.simplify.(J)

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

function generate_base_expressions_analytical(model::ContactModel)

	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c

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
	J = reshape(J, size(J_func(model, zeros(nq))))
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
	generate_dynamics_expressions(model::ContactModel)
Generate fast dynamics methods using Symbolics symbolic computing tools.
"""
function generate_dynamics_expressions(model::ContactModel; derivs = false)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nc = model.dim.c
	ncf = size(J_func(model, zeros(nq)))[1]

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

"""
	instantiate_base!(model,
		path::AbstractString="model/base.jld2")
Loads the base expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_base!(model::ContactModel, path::AbstractString="model/base.jld2")
	instantiate_base!(model.base, load_expressions(path))
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
function instantiate_dynamics!(model::ContactModel, path::AbstractString="model/dynamics.jld2";
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
