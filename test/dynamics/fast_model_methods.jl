@testset "Model Methods: Fast Base" begin
    dyn_path = joinpath(@__DIR__, "../../src/dynamics")
    include(joinpath(dyn_path, "particle/model.jl"))
	model = particle
    ContactControl.instantiate_base!(model, joinpath(dyn_path, "particle/base.jld2"))

    # Setup variables
    T = Float64
    nq = model.dim.q
    q0s = rand(SizedVector{nq,T})
    q̇0s = rand(SizedVector{nq,T})

    # Testing fast methods
    @test norm(ContactControl.M_fast(model, q0s) - M_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactControl.B_fast(model, q0s) - B_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactControl.N_fast(model, q0s) - N_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactControl.P_fast(model, q0s) - P_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactControl.C_fast(model, q0s, q̇0s) - C_func(model, q0s, q̇0s), Inf) < 1.0e-8
end

@testset "Model Methods: Fast Dynamics" begin
	dyn_path = joinpath(@__DIR__, "../../src/dynamics")
	include(joinpath(dyn_path, "particle/model.jl"))
	model = particle
	ContactControl.instantiate_dynamics!(model, joinpath(dyn_path, "particle/dynamics.jld2"))

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
    @test norm(ContactControl.d_fast(model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s) -
        dynamics(model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s), Inf) < 1.0e-8

    ContactControl.dq0_fast!(∇q0s, model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fq0(x) = dynamics(model, h, x, q1s, u1s, w1s, γ1s, b1s, q2s)
    @test norm(∇q0s - ForwardDiff.jacobian(fq0, q0s), Inf) < 1.0e-8

    ContactControl.dq1_fast!(∇q1s, model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fq1(x) = dynamics(model, h, q0s, x, u1s, w1s, γ1s, b1s, q2s)
    @test norm(∇q1s - ForwardDiff.jacobian(fq1, q1s), Inf) < 1.0e-8

    ContactControl.du1_fast!(∇u1s, model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fu1(x) = dynamics(model, h, q0s, q1s, x, w1s, γ1s, b1s, q2s)
    @test norm(∇u1s - ForwardDiff.jacobian(fu1, u1s), Inf) < 1.0e-8

    ContactControl.dγ1_fast!(∇γ1s, model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fγ1(x) = dynamics(model, h, q0s, q1s, u1s, w1s, x, b1s, q2s)
    @test norm(∇γ1s - ForwardDiff.jacobian(fγ1, γ1s), Inf) < 1.0e-8

    ContactControl.db1_fast!(∇b1s, model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fb1(x) = dynamics(model, h, q0s, q1s, u1s, w1s, γ1s, x, q2s)
    @test norm(∇b1s - ForwardDiff.jacobian(fb1, b1s), Inf) < 1.0e-8

    ContactControl.dq2_fast!(∇q2s, model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	fq2(x) = dynamics(model, h, q0s, q1s, u1s, w1s, γ1s, b1s, x)
    @test norm(∇q2s - ForwardDiff.jacobian(fq2, q2s), Inf) < 1.0e-8

    # dy_fast!(∇ys, model, h, q0s, q1s, u1s, w1s, γ1s, b1s, q2s)
	# # @show [∇q0s ∇q1s ∇u1s ∇γ1s ∇b1s ∇q2s]
    # @test norm(∇ys - [∇q0s ∇q1s ∇u1s ∇γ1s ∇b1s ∇q2s], Inf) < 1.0e-8
end

@testset "Model Methods: Fast Residual" begin
	res_path = joinpath(@__DIR__, "../../src/dynamics")
	include(joinpath(res_path, "particle/model.jl"))
	model = particle
	ContactControl.instantiate_base!(model, joinpath(res_path, "particle/base.jld2"))
	ContactControl.instantiate_dynamics!(model, joinpath(res_path, "particle/dynamics.jld2"))
	ContactControl.instantiate_residual!(model, joinpath(res_path, "particle/residual.jld2"))
	@load joinpath(res_path, "particle/sparse_jacobians.jld2") rz_sp rθ_sp

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
    ContactControl.r_fast!(rs, model, zs, θs, κs)
    @test norm(rs - residual(model, zs, θs, κs), Inf) < 1.0e-8

    ContactControl.rz_fast!(rz_sp, model, zs, θs, κs)
	fz(x) = residual(model, x, θs, κs)
    @test norm(rz_sp - ForwardDiff.jacobian(fz, zs), Inf) < 1.0e-8

	ContactControl.rθ_fast!(rθ_sp, model, zs, θs, κs)
	fθ(x) = residual(model, zs, x, κs)
    @test norm(rθ_sp - ForwardDiff.jacobian(fθ, θs), Inf) < 1.0e-8
end
