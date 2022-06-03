"""
	generate_base_expressions(model::Model)
Generate fast base methods using Symbolics symbolic computing tools.
"""
function generate_base_expressions(model::Model;
	M_analytical = true,
	C_analytical = true,
	mapping = (a, x) -> I,
	nv = model.nq)

	nq = model.nq
	nu = model.nu
	nw = model.nw
	nc = model.nc

	# Declare variables
	@variables q[1:nq]
	@variables q̇[1:nv]

	# Lagrangian
	L = lagrangian(model, q, q̇)
	L = Symbolics.simplify.(L)

	if !(M_analytical && C_analytical)
		ddL = Symbolics.hessian(L, [q; q̇], simplify=true)
	end

	# Mass Matrix
	if M_analytical
		M = M_func(model, q)
		M = reshape(M, (nv, nv))
		M = Symbolics.simplify.(M)
	else
		# ddL = Symbolics.hessian(L, [q; q̇], simplify=true)
		M = ddL[nq .+ (1:nv), nq .+ (1:nv)]
	end

	# Bias terms
	if C_analytical
		C = C_func(model, q, q̇)
		C = simplify.(C)
	else
		dLq = mapping(model, q)' * Symbolics.gradient(L, q, simplify=true) # including mapping for orientation (e.g., attitude Jacobian)
		# dLq̇ = Symbolics.gradient(L, q̇, simplify=true)
		# ddL = Symbolics.hessian(L, [q; q̇], simplify=true)
		ddLq̇q = ddL[nq .+ (1:nv), 1:nq] * mapping(model, q)

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

	# Kinematics
	k = kinematics(model, q)
	k = simplify.(k)

	# Build function
	expr = Dict{Symbol, Expr}()
	expr[:L]    = build_function([L], q, q̇, checkbounds=true)[1] # need to transpose to get a line vector
	expr[:M]    = build_function(M, q, checkbounds=true)[1]
	expr[:B]    = build_function(B, q, checkbounds=true)[1]
	expr[:A]    = build_function(A, q, checkbounds=true)[1]
	expr[:C]    = build_function(C, q, q̇, checkbounds=true)[1]
	expr[:k]    = build_function(k, q, checkbounds=true)[1]

	return expr
end

"""
	instantiate_base!(model,
		path::AbstractString="model/base.jld2")
Loads the base expressions from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_base!(model::Model, path::AbstractString="model/base.jld2")
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
	fct.C    = eval(expr[:C])
	fct.k    = eval(expr[:k])
	return nothing
end


"""
	generate_dynamics_expressions(model::Model)
Generate fast dynamics methods using Symbolics symbolic computing tools.
"""
function generate_dynamics_expressions(model::Model; derivs = false, nv = model.nq)
	nq = model.nq
	nu = model.nu
	nw = model.nw
	nc = model.nc
	# ncf = nc * dim(env)

	# Declare variables
	@variables q0[1:nq]
	@variables q1[1:nq]
	@variables u1[1:nu]
	@variables w1[1:nw]
	@variables Λ1[1:nv]
	@variables q2[1:nq]
	@variables h[1:1]

	# Expressions
	expr = Dict{Symbol, Expr}()

	# Dynamics
	d = dynamics(model, h, q0, q1, u1, w1, Λ1, q2)
	d = Symbolics.simplify.(d)

	# Functions
	expr[:d] = build_function(d, h, q0, q1, u1, w1, Λ1, q2, checkbounds=true)[1]

	if derivs
		∂q2 = Symbolics.jacobian(d, q2, simplify = true)
		dΛ1 = Symbolics.jacobian(d, Λ1, simplify = true)
		dθ =  Symbolics.jacobian(d, [q0; q1; u1; w1; h])#, simplify = true)

		expr[:∂q2] = build_function(∂q2, h, q0, q1, u1, w1, Λ1, q2, checkbounds=true)[1]
		expr[:dΛ1] = build_function(dΛ1, h, q0, q1, u1, w1, Λ1, q2, checkbounds=true)[1]
		expr[:dθ]  = build_function(dθ,  h, q0, q1, u1, w1, Λ1, q2, checkbounds=true)[1]
	end

	return expr
end

"""
	instantiate_dynamics!(model,
		path::AbstractString="dynamics.jld2")
Loads the dynamics expressoins from the `path`, evaluates them to generate functions,
stores them into the model.
"""
function instantiate_dynamics!(model::Model, path::AbstractString="model/dynamics.jld2";
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
		fct.dΛ1 = eval(expr[:dΛ1])
		fct.∂q2 = eval(expr[:∂q2])
	end

	return nothing
end
