@testset "Fast Base Model Methods" begin
    dyn_path = joinpath(@__DIR__, "../../src/dynamics")
    include(joinpath(dyn_path, "quadruped_model.jl"))
    model = deepcopy(quadruped)
    # expr_bas = generate_base_expressions(model)
    # save_expressions(expr_bas, joinpath(dyn_path, ".expr/quadruped_base.jld2"), overwrite=true)
    instantiate_base!(model, joinpath(dyn_path, ".expr/quadruped_base.jld2"))

    # Setup variables
    T = Float64
    nq = model.dim.q
    q0s = rand(SizedVector{nq,T})
    q̇0s = rand(SizedVector{nq,T})

    # Testing fast methods
    @test norm(M_fast(model, q0s) - M_func(model, q0s), 1) < 1e-14
    @test norm(B_fast(model, q0s) - B_func(model, q0s), 1) < 1e-14
    @test norm(N_fast(model, q0s) - N_func(model, q0s), 1) < 1e-14
    @test norm(P_fast(model, q0s) - P_func(model, q0s), 1) < 1e-14
    @test norm(C_fast(model, q0s, q̇0s) - _C_func(model, q0s, q̇0s), 1) < 1e-12
end

@testset "Fast Dynamics Model Methods" begin
    dyn_path = joinpath(@__DIR__, "../../src/dynamics")
    include(joinpath(dyn_path, "quadruped_model.jl"))
    model = deepcopy(quadruped)
    instantiate_base!(model, joinpath(dyn_path, ".expr/quadruped_base.jld2"))
    # expr_dyn = generate_dynamics_expressions(model)
    # save_expressions(expr_dyn, joinpath(dyn_path, ".expr/quadruped_dynamics.jld2"), overwrite=true)
    instantiate_dynamics!(model, joinpath(dyn_path, ".expr/quadruped_dynamics.jld2"))

    # Setup variables
    T = Float64
    nq = model.dim.q
    nu = model.dim.u
    nγ = model.dim.γ
    nb = model.dim.b
    ny = model.dim.y

    q0s = rand(SizedVector{nq,T})
    q1s = rand(SizedVector{nq,T})
    u1s = rand(SizedVector{nu,T})
    γ1s = rand(SizedVector{nγ,T})
    b1s = rand(SizedVector{nb,T})
    q2s = rand(SizedVector{nq,T})

    ∇ys   = rand(SizedMatrix{nq,ny,T})
    ∇q0s  = rand(SizedMatrix{nq,nq,T})
    ∇q1s  = rand(SizedMatrix{nq,nq,T})
    ∇u1s  = rand(SizedMatrix{nq,nu,T})
    ∇γ1s  = rand(SizedMatrix{nq,nγ,T})
    ∇b1s  = rand(SizedMatrix{nq,nb,T})
    ∇q2s  = rand(SizedMatrix{nq,nq,T})

    # Testing dynamics methods
    @test norm(dynamics_fast(model, q0s, q1s, u1s, γ1s, b1s, q2s) -
        dynamics(model, model.dt, q0s, q1s, u1s, γ1s, b1s, q2s), 1) < 1e-12

    ∇q0_dynamics_fast!(model, ∇q0s, q0s, q1s, u1s, γ1s, b1s, q2s)
    @test norm(∇q0s - ∇q0_dynamics(model, model.dt, q0s, q1s, u1s, γ1s, b1s, q2s), 1) < 1e-12

    ∇q1_dynamics_fast!(model, ∇q1s, q0s, q1s, u1s, γ1s, b1s, q2s)
    @test norm(∇q1s - ∇q1_dynamics(model, model.dt, q0s, q1s, u1s, γ1s, b1s, q2s), 1) < 1e-12

    ∇u1_dynamics_fast!(model, ∇u1s, q0s, q1s, u1s, γ1s, b1s, q2s)
    @test norm(∇u1s - ∇u1_dynamics(model, model.dt, q0s, q1s, u1s, γ1s, b1s, q2s), 1) < 1e-12

    ∇γ1_dynamics_fast!(model, ∇γ1s, q0s, q1s, u1s, γ1s, b1s, q2s)
    @test norm(∇γ1s - ∇γ1_dynamics(model, model.dt, q0s, q1s, u1s, γ1s, b1s, q2s), 1) < 1e-12

    ∇b1_dynamics_fast!(model, ∇b1s, q0s, q1s, u1s, γ1s, b1s, q2s)
    @test norm(∇b1s - ∇b1_dynamics(model, model.dt, q0s, q1s, u1s, γ1s, b1s, q2s), 1) < 1e-12

    ∇q2_dynamics_fast!(model, ∇q2s, q0s, q1s, u1s, γ1s, b1s, q2s)
    @test norm(∇q2s - ∇q2_dynamics(model, model.dt, q0s, q1s, u1s, γ1s, b1s, q2s), 1) < 1e-12

    ∇y_dynamics_fast!(model, ∇ys, q0s, q1s, u1s, γ1s, b1s, q2s)
    @test norm(∇ys - [∇q0s ∇q1s ∇u1s ∇γ1s ∇b1s ∇q2s], 1) < 1e-12
end

@testset "Fast Residual Model Methods" begin
    dyn_path = joinpath(@__DIR__, "../../src/dynamics")
    include(joinpath(dyn_path, "quadruped_model.jl"))
    model = deepcopy(quadruped)
    instantiate_base!(model, joinpath(dyn_path, ".expr/quadruped_base.jld2"))
    instantiate_dynamics!(model, joinpath(dyn_path, ".expr/quadruped_dynamics.jld2"))
    # expr_res = generate_residual_expressions(model)
    # save_expressions(expr_dyn, joinpath(dyn_path, ".expr/quadruped_residual.jld2"), overwrite=true)
    instantiate_residual!(model, joinpath(dyn_path, ".expr/quadruped_residual.jld2"))

    # Setup variables
    T = Float64
    nz = model.dim.z
    nθ = model.dim.θ

    zs = rand(SizedVector{nz})
    θs = rand(SizedVector{nθ})
    ρs = 1e-3
    rs = rand(SizedVector{nz})
    ∇zs = rand(nz,nz)
    ∇θs = rand(nz,nθ)

    # Testing residual methods
    r_fast!(model, rs, zs, θs, ρs)
    @test norm(rs - residual(model, model.dt, zs, θs, ρs), 1) < 1e-12

    rz_fast!(model, ∇zs, zs, θs, ρs)
    @test norm(∇zs - ∇z_residual(model, model.dt, zs, θs, ρs), 1) < 1e-12

    rθ_fast!(model, ∇θs, zs, θs, ρs)
    @test norm(∇θs - ∇θ_residual(model, model.dt, zs, θs, ρs), 1) < 1e-12
end
