@testset "Solver: Random QP w/ equalities (QDLDL)" begin
    """
        minimize   x' P x + q' x
        subject to    Ax = b
    """

    #####
    n = 10
    m = 3

    nw = n + m
    nθ = 2n + m * n + m

    function obj(z, θ)
        x = z[1:n]

        P = Diagonal(θ[1:n])
        p = θ[n .+ (1:n)]
        return transpose(x) * P * x + transpose(p) * x
    end

    function constraints(z, θ)
        x = z[1:n]
        A = reshape(θ[2n .+ (1:(m * n))], m, n)
        b = θ[2n + m * n .+ (1:m)]

        A * x - b
    end

    function lagrangian(w, θ)
        z = w[1:n]
        y = w[n .+ (1:m)]

        L = 0.0
        L += obj(z, θ)
        L += dot(y, constraints(z, θ))

        return L
    end

    # residual
    function _r(w, θ, κ)
        L = lagrangian(w, θ)
        Lz = Symbolics.gradient(L, w[1:n])
        c = constraints(w[1:n], θ)

        [
            Lz;
            c;
        ]
    end

    @variables z_sym[1:nw]
    @variables θ_sym[1:nθ]
    @variables κ_sym[1:1]

    r_sym = _r(z_sym, θ_sym, κ_sym)
    rf! = eval(Symbolics.build_function(r_sym, z_sym, θ_sym, κ_sym)[2])
    rz_exp = Symbolics.jacobian(r_sym, z_sym)
    rθ_exp = Symbolics.jacobian(r_sym, θ_sym)
    rz_sp = similar(rz_exp, Float64)
    rθ_sp = similar(rθ_exp, Float64)
    rzf! = eval(Symbolics.build_function(rz_exp, z_sym, θ_sym)[2])
    rθf! = eval(Symbolics.build_function(rθ_exp, z_sym, θ_sym)[2])

    x0 = randn(n)
    A = rand(m, n)
    b = A * x0
    z = [randn(n); zeros(m)]
    θ = [ones(n); zeros(n); vec(A); b]
    κ = [1.0]

    r = zeros(nw)
    rz = zeros(nw, nw)
    rθ = zeros(nw, nθ)

    opts = ContactImplicitMPC.InteriorPointOptions(diff_sol=true,
        undercut=10.0,
        max_ls=25,
        max_iter=100,
        verbose=false,
        solver=:lu_solver)

    idx = ContactImplicitMPC.IndicesOptimization(
        nw, nw,
        [collect(1:0), collect(1:0)],
        [collect(1:0), collect(1:0)],
        Vector{Vector{Int}}[],
        Vector{Vector{Int}}[],
        collect(1:(n + m)),
        Vector{Int}(),
        Vector{Int}(),
        Vector{Int}[],
        Vector{Int}())

    # solver
    ip = ContactImplicitMPC.interior_point(z, θ,
        idx=idx,
        r! = rf!, rz! = rzf!, rθ! = rθf!,
        rz = rz,
        rθ = rθ,
        opts = opts)

    status = ContactImplicitMPC.interior_point_solve!(ip)
    # @benchmark CALIPSO.interior_point_solve!($ip)

    # test
    @test status
    @test norm(ip.r, Inf) < opts.r_tol
    @test norm(A * ip.z[1:n] - b, Inf) < 1.0e-6
    @test norm(ip.δz, 1) != 0.0
end
