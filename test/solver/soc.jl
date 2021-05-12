@testset "Solver: Second-order cone program (friction)" begin
	"""
		minimize   v' b
		subject to ||b|| <= μn
	"""

	n = 7
	m = 3

	idx_ineq = collect(1:0)
	idx_soc = [collect(1:3), collect(5:7)]

	# residual
	function _r!(r, z, θ, κ)
		x = z[1:3]
		y = z[4:4]
		z = z[5:7]

		v = θ[1:2]
		μγ = θ[3]

		r[1:3] = [0.0; v] + [y; 0.0; 0.0] - z
		r[4] = x[1] - μγ
		r[5:7] = second_order_cone_product(x, z) - κ * [1.0; 0.0; 0.0]
		nothing
	end

	@variables r_sym[1:n]
	@variables z_sym[1:n]
	@variables θ_sym[1:m]
	@variables κ_sym

	parallel = Symbolics.SerialForm()
	_r!(r_sym, z_sym, θ_sym, κ_sym)
	r_sym = simplify.(r_sym)
	rf! = eval(Symbolics.build_function(r_sym, z_sym, θ_sym, κ_sym,
		parallel = parallel)[2])
	rz_exp = Symbolics.jacobian(r_sym, z_sym, simplify = true)
	rθ_exp = Symbolics.jacobian(r_sym, θ_sym, simplify = true)
	rz_sp = similar(rz_exp, Float64)
	rθ_sp = similar(rθ_exp, Float64)
	rzf! = eval(Symbolics.build_function(rz_exp, z_sym, θ_sym,
		parallel = parallel)[2])
	rθf! = eval(Symbolics.build_function(rθ_exp, z_sym, θ_sym,
		parallel = parallel)[2])

	# problem setup
	v = [0.0; 0.0]
	μγ = 1.0
	θ = [v; μγ]

	z = ones(n)
	z[1] += 1.0
	z[5] += 1.0

	# solver
	ip = ContactControl.interior_point(z, θ,
		idx_ineq = idx_ineq,
		idx_soc = idx_soc,
		r! = rf!, rz! = rzf!, rθ! = rθf!,
		rz = rz_sp,
		rθ = rθ_sp,
		opts = ContactControl.InteriorPointOptions(diff_sol = true))

	# solve
	status = ContactControl.interior_point!(ip)

	# test
	@test status
	@test norm(ip.r, Inf) < opts.r_tol
	@test ContactControl.inequality_check(ip.z, ip.idx_ineq)
	@test ip.κ[1] < opts.κ_tol
	@test norm(ip.δz, 1) != 0.0

	# problem setup
	v = [1.0; 0.0]
	μγ = 1.0
	θ = [v; μγ]

	z = ones(n)
	z[1] += 1.0
	z[5] += 1.0

	ip.z .= z
	ip.θ .= θ

	# solve
	status = ContactControl.interior_point!(ip)

	# test
	@test status
	@test norm(ip.r, Inf) < opts.r_tol
	@test ContactControl.inequality_check(ip.z, ip.idx_ineq)
	@test ip.κ[1] < opts.κ_tol
	@test norm(ip.δz, 1) != 0.0

	# problem setup
	v = [1.0; 10.0]
	μγ = 1.0
	θ = [v; μγ]

	z = ones(n)
	z[1] += 1.0
	z[5] += 1.0

	ip.z .= z
	ip.θ .= θ

	# solve
	status = ContactControl.interior_point!(ip)

	# test
	@test status
	@test norm(ip.r, Inf) < opts.r_tol
	@test ContactControl.inequality_check(ip.z, ip.idx_ineq)
	@test ip.κ[1] < opts.κ_tol
	@test norm(ip.δz, 1) != 0.0

	# problem setup
	v = [1.0; 10.0]
	μγ = 0.0
	θ = [v; μγ]

	z = ones(n)
	z[1] += 1.0
	z[5] += 1.0

	ip.z .= z
	ip.θ .= θ

	# solve
	status = ContactControl.interior_point!(ip)

	# test
	@test status
	@test norm(ip.r, Inf) < opts.r_tol
	@test ContactControl.inequality_check(ip.z, ip.idx_ineq)
	@test ip.κ[1] < opts.κ_tol
	@test norm(ip.δz, 1) != 0.0

	# problem setup
	v = [0.0; 0.0]
	μγ = 0.0
	θ = [v; μγ]

	z = ones(n)
	z[1] += 1.0
	z[5] += 1.0

	ip.z .= z
	ip.θ .= θ

	# solve
	status = ContactControl.interior_point!(ip)

	# test
	@test status
	@test norm(ip.r, Inf) < opts.r_tol
	@test ContactControl.inequality_check(ip.z, ip.idx_ineq)
	@test ip.κ[1] < opts.κ_tol
	@test norm(ip.δz, 1) != 0.0
end
