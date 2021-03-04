@testset "Fast Model Methods" begin
    dyn_path = joinpath(@__DIR__, "../../src/dynamics")
    include(joinpath(dyn_path, "quadruped_model.jl"))
    model = deepcopy(quadruped)
    # expr_bas = generate_base_expressions(model)
    # save_expressions(expr_bas, joinpath(@__DIR__, ".expr/quadruped_base.jld2"), overwrite=true)
    instantiate_base!(model, joinpath(dyn_path, ".expr/quadruped_base.jld2"))

    T = Float64
    nq = model.dim.q
    q0s = rand(SizedVector{nq,T})
    q̇0s = rand(SizedVector{nq,T})

    # Testing fast methods
    @test norm(M_fast(model, q0s) - M_func(model, q0s), 1) < 1e-14
    @test norm(B_fast(model, q0s) - B_func(model, q0s), 1) < 1e-14
    @test norm(N_fast(model, q0s) - N_func(model, q0s), 1) < 1e-14
    @test norm(P_fast(model, q0s) - _P_func(model, q0s), 1) < 1e-14
    @test norm(C_fast(model, q0s, q̇0s) - _C_func(model, q0s, q̇0s), 1) < 1e-12
end
