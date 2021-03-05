@testset "Fast Base Model Methods" begin
    include(joinpath(pwd(), "src/dynamics/quadruped/model.jl"))

    # Setup variables
    T = Float64
    nq = model.dim.q
    q0s = rand(SizedVector{nq,T})
    q̇0s = rand(SizedVector{nq,T})

    # Testing fast methods
    @test norm(M_fast(q0s) - M_func(model, q0s), 1) < 1.0e-14
    @test norm(B_fast(q0s) - B_func(model, q0s), 1) < 1.0e-14
    @test norm(N_fast(q0s) - N_func(model, q0s), 1) < 1.0e-14
    @test norm(P_fast(q0s) - P_func(model, q0s), 1) < 1.0e-14
    @test norm(C_fast(q0s, q̇0s) - C_func(model, q0s, q̇0s), 1) < 1.0e-12
end

@testset "Fast Dynamics Model Methods" begin
    include(joinpath(pwd(), "src/dynamics/quadruped/model.jl"))

    # Setup variables
    T = Float64
	h = 0.1
    nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
    nb = model.dim.b

    q0s = rand(SizedVector{nq,T})
    q1s = rand(SizedVector{nq,T})
    u1s = rand(SizedVector{nu,T})
    u1s = rand(SizedVector{nu,T})
    w1s = rand(SizedVector{nw,T})
    γ1s = rand(SizedVector{nc,T})
    b1s = rand(SizedVector{nb,T})
    q2s = rand(SizedVector{nq,T})

    ∇ys   = rand(SizedMatrix{nq,2nq + nu + nc + nb + nq,T})
    ∇q0s  = rand(SizedMatrix{nq,nq,T})
    ∇q1s  = rand(SizedMatrix{nq,nq,T})
    ∇u1s  = rand(SizedMatrix{nq,nu,T})
    ∇γ1s  = rand(SizedMatrix{nq,nc,T})
    ∇b1s  = rand(SizedMatrix{nq,nb,T})
    ∇q2s  = rand(SizedMatrix{nq,nq,T})

    # Testing dynamics methods
    @test norm(d(h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s) -
        dynamics(model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s), 1) < 1e-12

    dq0!(∇q0s, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fq0(x) = dynamics(model, h, x, q1s, u1s, w1s, γ1s, b1s, q2s)
    @test norm(∇q0s - ForwardDiff.jacobian(fq0, q0s), 1) < 1e-12

    dq1!(∇q1s, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fq1(x) = dynamics(model, h, q0s, x, u1s, w1s, γ1s, b1s, q2s)
    @test norm(∇q1s - ForwardDiff.jacobian(fq1, q1s), 1) < 1e-12

    du1!(∇u1s, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fu1(x) = dynamics(model, h, q0s, q1s, x, w1s, γ1s, b1s, q2s)
    @test norm(∇u1s - ForwardDiff.jacobian(fu1, u1s), 1) < 1e-12

    dγ1!(∇γ1s, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fγ1(x) = dynamics(model, h, q0s, q1s, u1s, w1s, x, b1s, q2s)
    @test norm(∇γ1s - ForwardDiff.jacobian(fγ1, γ1s), 1) < 1e-12

    db1!(∇b1s, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fb1(x) = dynamics(model, h, q0s, q1s, u1s, w1s, γ1s, x, q2s)
    @test norm(∇b1s - ForwardDiff.jacobian(fb1, b1s), 1) < 1e-12

    dq2!(∇q2s, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fq2(x) = dynamics(model, h, q0s, q1s, u1s, w1s, γ1s, b1s, x)
    @test norm(∇q2s - ForwardDiff.jacobian(fq2, q2s), 1) < 1e-12

    dy!(∇ys, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
    @test norm(∇ys - [∇q0s ∇q1s ∇u1s ∇γ1s ∇b1s ∇q2s], 1) < 1e-12
end

@testset "Fast Residual Model Methods" begin
    include(joinpath(pwd(), "src/dynamics/quadruped/model.jl"))

    # Setup variables
    T = Float64
	h = 0.1
    nz = num_var(model)
    nθ = num_data(model)

    zs = rand(SizedVector{nz})
    θs = rand(SizedVector{nθ})
    κs = 1.0
    rs = rand(SizedVector{nz})

    # Testing residual methods
    r!(rs, zs, θs, κs)
    @test norm(rs - residual(model, zs, θs, κs), 1) < 1e-12

    rz!(rz_sp, zs, θs, κs)
	fz(x) = residual(model, x, θs, κs)
    @test norm(rz_sp - ForwardDiff.jacobian(fz, zs), 1) < 1e-12

	rθ!(rθ_sp, zs, θs, κs)
	fθ(x) = residual(model, zs, x, κs)
    @test norm(rθ_sp - ForwardDiff.jacobian(fθ, θs), 1) < 1e-12
end
