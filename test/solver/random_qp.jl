@testset "Solver: Random QP (LU)" begin
    """
        minimize   x' P x + q' x
        subject to    x >= 0
    """

    n = 1000

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
        x = z[1:1000]
        y = z[1001:2000]
        P = θ[1:1000]
        q = θ[1001:2000]

        r[1:1000] = 2.0 * Diagonal(P) * x + q - y
        r[1001:2000] = Diagonal(x) * y .- κ
        nothing
    end

    @variables r_sym[1:2000]
    @variables z_sym[1:2000]
    @variables θ_sym[1:2000]
    @variables κ_sym[1:1]

    parallel = false
    _r!(r_sym, z_sym, θ_sym, κ_sym)
    r_sym = simplify.(r_sym)
    rf! = eval(ModelingToolkit.build_function(r_sym, z_sym, θ_sym, κ_sym,
        parallel = parallel)[2])
    rz_exp = ModelingToolkit.sparsejacobian(r_sym, z_sym, simplify = true)
    rθ_exp = ModelingToolkit.sparsejacobian(r_sym, θ_sym, simplify = true)
    rz_sp = similar(rz_exp, Float64)
    rθ_sp = similar(rθ_exp, Float64)
    rzf! = eval(ModelingToolkit.build_function(rz_exp, z_sym, θ_sym,
        parallel = parallel)[2])
    rθf! = eval(ModelingToolkit.build_function(rθ_exp, z_sym, θ_sym,
        parallel = parallel)[2])

    # solver
    ip = ContactControl.interior_point(2 * n, 2 * n, idx_ineq,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp)

    # options
    opts = ContactControl.InteriorPointOptions(diff_sol = true)

    # solve
    status = ContactControl.interior_point!(ip, z, θ, opts = opts)

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
    @variables κ_sym[1:1]

    parallel = false
    _r!(r_sym, z_sym, θ_sym, κ_sym)
    r_sym = simplify.(r_sym)
    rf! = eval(ModelingToolkit.build_function(r_sym, z_sym, θ_sym, κ_sym,
        parallel = parallel)[2])
    rz_exp = ModelingToolkit.sparsejacobian(r_sym, z_sym, simplify = true)
    # rθ_exp = ModelingToolkit.jacobian(r_sym, θ_sym, simplify = false)
    rz_sp = similar(rz_exp, Float64)
    rθ_sp = zeros(0, 0) #similar(rθ_exp, Float64)
    rzf! = eval(ModelingToolkit.build_function(rz_exp, z_sym, θ_sym,
        parallel = parallel)[2])
    rθf! = x -> nothing #eval(ModelingToolkit.build_function(rθ_exp, z_sym, θ_sym,
        # parallel = parallel)[2])

    # solver
    ip = ContactControl.interior_point(n + m, 2 * n + m * n, idx_ineq,
        x = z, θ = θ,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz_sp,
        rθ = rθ_sp,
        solver = :ldl_solver)

    opts = ContactControl.InteriorPointOptions(diff_sol = false,
        κ_init = 1.0,
        κ_tol = 1.0,
        max_iter_outer = 1)

    status = ContactControl.interior_point!(ip, z, θ, opts = opts)

    # test
    @test status
    @test norm(ip.r, Inf) < opts.r_tol
    @test norm(A * ip.z[1:n] - b, Inf) < opts.r_tol
end
