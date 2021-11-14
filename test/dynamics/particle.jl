@testset "Model Methods: Fast Base (Particle)" begin
    s = ContactImplicitMPC.get_simulation("particle", "flat_3D_lc", "flat_lc")
	model = s.model

    # Setup variables
    T = Float64
    nq = model.nq
    q0s = rand(nq)
    q̇0s = rand(nq) # TODO: Symbolics not happy with SizedArray

    # Testing fast methods
    @test norm(ContactImplicitMPC.M_fast(model, q0s) - ContactImplicitMPC.M_func(model, q0s), Inf) < 1.0e-8
	@test norm(ContactImplicitMPC.B_fast(model, q0s) - ContactImplicitMPC.B_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactImplicitMPC.A_fast(model, q0s) - ContactImplicitMPC.A_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactImplicitMPC.C_fast(model, q0s, q̇0s) - ContactImplicitMPC.C_func(model, q0s, q̇0s), Inf) < 1.0e-8
end

@testset "Model Methods: Fast Dynamics (Particle)" begin
	s = ContactImplicitMPC.get_simulation("particle", "flat_3D_lc", "flat_lc")
    model = s.model
	env = s.env

	# Setup variables
    T = Float64
    nq = model.nq
    nu = model.nu
    nw = model.nw
    nc = model.nc
	nb = nc * friction_dim(s.env)

	hs = 0.1 * ones(1)
    q0s = rand(nq)
    q1s = rand(nq)
    u1s = rand(nu)
    u1s = rand(nu)
    w1s = rand(nw)
    γ1s = rand(nc)
    b1s = rand(nb)
	λ1s = rand(nc * dim(env))
    q2s = rand(nq)

    # Testing dynamics methods
    @test norm(ContactImplicitMPC.d_fast(model, hs, q0s, q1s, u1s, w1s, λ1s, q2s) -
        ContactImplicitMPC.dynamics(model, hs, q0s, q1s, u1s, w1s, λ1s, q2s), Inf) < 1.0e-8
end
