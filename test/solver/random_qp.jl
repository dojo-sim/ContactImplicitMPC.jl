@testset "Solver: Random QP (LU)" begin
    """
        minimize   x' P x + q' x
        subject to    x >= 0
    """

    n = 100

    _P = rand(n)
    P = Diagonal(_P)
    q = rand(n)
    θ = [_P; q]
    z = ones(2 * n)
    r = zeros(2 * n)
    rz = zeros(2 * n, 2 * n)
    rθ = zeros(2 * n, 2 * n)
    κ = 1.0

    idx_ineq = collect(1:2 * n)

    # residual
    function _r!(r, z, θ, κ)
        x = z[1:100]
        y = z[101:200]
        P = θ[1:100]
        q = θ[101:200]

        r[1:100] = 2.0 * Diagonal(P) * x + q - y
        r[101:200] = Diagonal(x) * y .- κ
        nothing
    end

    @variables r_sym[1:200]
    @variables z_sym[1:200]
    @variables θ_sym[1:200]
    @variables κ_sym

    parallel = Symbolics.SerialForm()
    _r!(r_sym, z_sym, θ_sym, κ_sym)
    r_sym = simplify.(r_sym)
    rf! = eval(Symbolics.build_function(r_sym, z_sym, θ_sym, κ_sym,
        parallel = parallel)[2])
    rz_exp = Symbolics.sparsejacobian(r_sym, z_sym, simplify = true)
    rθ_exp = Symbolics.sparsejacobian(r_sym, θ_sym, simplify = true)
    rz_sp = similar(rz_exp, Float64)
    rθ_sp = similar(rθ_exp, Float64)
    rzf! = eval(Symbolics.build_function(rz_exp, z_sym, θ_sym,
        parallel = parallel)[2])
    rθf! = eval(Symbolics.build_function(rθ_exp, z_sym, θ_sym,
        parallel = parallel)[2])

    # options
    opts = ContactControl.InteriorPointOptions(diff_sol = true)

    # solver
    ip = ContactControl.interior_point(z, θ,
        idx_ineq = idx_ineq,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp,
        opts = opts)

    # solve
    status = ContactControl.interior_point!(ip)

    # test
    @test status
    @test norm(ip.r, Inf) < opts.r_tol
    @test !ContactControl.inequality_check(ip.z, ip.idx_ineq)
    @test ip.κ[1] < opts.κ_tol
    @test norm(ip.δz, 1) != 0.0
end

@testset "Solver: Random QP (LDL)" begin
    """
        minimize   x' P x + q' x
        subject to    A x = b
    """

    n = 10
    m = 5
    x = rand(n)
    A = randn(m, n)
    A_vec = vec(A)
    b = A * x
    _P = rand(n)
    P = Diagonal(_P)
    q = rand(n)
    θ = [_P; q; A_vec]
    z = ones(n + m)
    r = zeros(n + m)
    rz = zeros(n + m, n + m)
    rθ = zeros(n + m, 2 * n + length(A_vec))
    κ = 1.0

    idx_ineq = collect(1:0)

    # residual
    function _r!(r, z, θ, κ)
        x = z[1:10]
        y = z[11:15]
        P = θ[1:10]
        q = θ[11:20]
        Ā = reshape(θ[21:70], (5, 10))

        r[1:10] = 2.0 * Diagonal(P) * x + q + Ā' * y
        r[11:15] = Ā * x - b
        nothing
    end

    @variables r_sym[1:15]
    @variables z_sym[1:15]
    @variables θ_sym[1:70]
    @variables κ_sym

    parallel = Symbolics.SerialForm()
    _r!(r_sym, z_sym, θ_sym, κ_sym)
    r_sym = simplify.(r_sym)
    rf! = eval(Symbolics.build_function(r_sym, z_sym, θ_sym, κ_sym,
        parallel = parallel)[2])
    rz_exp = Symbolics.sparsejacobian(r_sym, z_sym, simplify = true)
    # rθ_exp = Symbolics.jacobian(r_sym, θ_sym, simplify = false)
    rz_sp = similar(rz_exp, Float64)
    rθ_sp = zeros(0, 0) #similar(rθ_exp, Float64)
    rzf! = eval(Symbolics.build_function(rz_exp, z_sym, θ_sym)[2])
    rθf! = x -> nothing #eval(Symbolics.build_function(rθ_exp, z_sym, θ_sym,
        # parallel = parallel)[2])

    # options
    opts = ContactControl.InteriorPointOptions(
        diff_sol = false,
        κ_init = 1.0,
        κ_tol = 1.0,
        max_iter_outer = 1,
        solver = :ldl_solver)

    # solver
    ip = ContactControl.interior_point(z, θ,
        idx_ineq = idx_ineq,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp,
        opts = opts)

    status = ContactControl.interior_point!(ip, z, θ)

    # test
    @test status
    @test norm(ip.r, Inf) < opts.r_tol
    @test norm(A * ip.z[1:n] - b, Inf) < opts.r_tol
end

@testset "Solver: Random QP (QR)" begin
    """
        minimize   x' P x + q' x
        subject to    x >= 0
    """

    n = 10

    _P = rand(n)
    P = Diagonal(_P)
    q = rand(n)
    θ = [_P; q]
    z = ones(2 * n)
    r = zeros(2 * n)
    rz = zeros(2 * n, 2 * n)
    rθ = zeros(2 * n, 2 * n)
    κ = 1.0

    idx_ineq = collect(1:2 * n)

    # residual
    function _r!(r, z, θ, κ)
        x = z[1:10]
        y = z[11:20]
        P = θ[1:10]
        q = θ[11:20]

        r[1:10] = 2.0 * Diagonal(P) * x + q - y
        r[11:20] = Diagonal(x) * y .- κ
        nothing
    end

    @variables r_sym[1:20]
    @variables z_sym[1:20]
    @variables θ_sym[1:20]
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

    # options
    opts_cgs = ContactControl.InteriorPointOptions(
        diff_sol = false, solver = :cgs_solver)

    opts_mgs = ContactControl.InteriorPointOptions(
        diff_sol = false, solver = :mgs_solver)

    opts_dmgs = ContactControl.InteriorPointOptions(
        diff_sol = false, solver = :dmgs_solver)

    # solver
    ip_cgs = ContactControl.interior_point(z, θ,
        idx_ineq = idx_ineq,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp,
        opts = opts_cgs)

    ip_mgs = ContactControl.interior_point(z, θ,
        idx_ineq = idx_ineq,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp,
        opts = opts_mgs)

    ip_dmgs = ContactControl.interior_point(z, θ,
        idx_ineq = idx_ineq,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp,
        opts = opts_dmgs)

    # solve
    status_cgs = ContactControl.interior_point!(ip_cgs, z, θ)
    status_mgs = ContactControl.interior_point!(ip_mgs, z, θ)
    status_dmgs = ContactControl.interior_point!(ip_dmgs, z, θ)

    # test
    @test status_cgs
    @test norm(ip_cgs.r, Inf) < opts_cgs.r_tol
    @test !ContactControl.inequality_check(ip_cgs.z, ip_cgs.idx_ineq)
    @test ip_cgs.κ[1] < opts_cgs.κ_tol

    @test status_mgs
    @test norm(ip_mgs.r, Inf) < opts_mgs.r_tol
    @test !ContactControl.inequality_check(ip_mgs.z, ip_mgs.idx_ineq)
    @test ip_mgs.κ[1] < opts_mgs.κ_tol

    @test status_dmgs
    @test norm(ip_dmgs.r, Inf) < opts_dmgs.r_tol
    @test !ContactControl.inequality_check(ip_dmgs.z, ip_dmgs.idx_ineq)
    @test ip_dmgs.κ[1] < opts_dmgs.κ_tol
end

# Compare with OSQP
# n = 10
# m = 5
#
# P = sparse(Diagonal(rand(n)))
# q = randn(n)
# x = randn(n)
# A = sparse(randn(m, n))
# b = A * x
# l = b - rand(m)
# u = b + rand(m)
#
# using OSQP
# model = OSQP.Model()
# options = Dict(:verbose => false,
#                :eps_abs => 1e-09,
#                :eps_rel => 1e-09,
#                :check_termination => 1,
#                :polish => true,
#                :max_iter => 4000,
#                :rho => 0.1,
#                :adaptive_rho => false,
#                :warm_start => true)
# OSQP.setup!(model, P = P, q = q, A = A, l = l, u = u)
#
# @time results = OSQP.solve!(model)
# u - A * results.x
# A * results.x - l
#
# # residual
# function _r!(r, z, θ, κ)
#     x = z[1:10]
#     s1 = z[11:15]
#     s2 = z[16:20]
#     y1 = z[21:25]
#     y2 = z[26:30]
#
#     r[1:10] = P * x + q + transpose(A) * y1 - transpose(A) * y2
#     r[11:15] = s1 - u + A * x
#     r[16:20] = s2 - A * x + l
#     r[21:25] = Diagonal(s1) * y1 .- κ
#     r[26:30] = Diagonal(s2) * y2 .- κ
#     nothing
# end
#
# @variables r_sym[1:30]
# @variables z_sym[1:30]
# @variables θ_sym
# @variables κ_sym
#
# parallel = Symbolics.SerialForm()
# _r!(r_sym, z_sym, θ_sym, κ_sym)
# r_sym = simplify.(r_sym)
# rf! = eval(Symbolics.build_function(r_sym, z_sym, θ_sym, κ_sym,
#     parallel = parallel)[2])
# rz_exp = Symbolics.jacobian(r_sym, z_sym, simplify = true)
# rθ_exp = x -> nothing #Symbolics.sparsejacobian(r_sym, θ_sym, simplify = true)
# rz_sp = similar(rz_exp, Float64)
# rθ_sp = zeros(0, 0) #similar(rθ_exp, Float64)
# rzf! = eval(Symbolics.build_function(rz_exp, z_sym, θ_sym,
#     parallel = parallel)[2])
# rθf! = x -> nothing #eval(Symbolics.build_function(rθ_exp, z_sym, θ_sym,
#     # parallel = parallel)[2])
#
# opts = ContactControl.InteriorPointOptions(
#     diff_sol = false, solver = :dmgs_solver)
#
# idx_ineq = collect(11:30)
#
# # solver
# z = ones(30)
# θ = zeros(0)
# ip = ContactControl.interior_point(z, θ,
#     idx_ineq = idx_ineq,
#     r! = rf!, rz! = rzf!, rθ! = rθf!,
#     rz = rz_sp,
#     rθ = rθ_sp,
#     opts = opts)
#
# linear_solve!(ip.solver, ip.Δ, ip.rz, ip.r)
# # solve
# @time status = ContactControl.interior_point!(ip, copy(z), copy(θ))
# # norm(ip.z[1:10] - results.x)
