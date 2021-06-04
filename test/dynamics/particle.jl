@testset "Model Methods: Fast Base (Particle)" begin
    s = ContactControl.get_simulation("particle", "flat_3D_lc", "flat")
	model = s.model

    # Setup variables
    T = Float64
    nq = model.dim.q
    q0s = rand(nq)
    q̇0s = rand(nq) # TODO: Symbolics not happy with SizedArray

    # Testing fast methods
    @test norm(ContactControl.M_fast(model, q0s) - ContactControl.M_func(model, q0s), Inf) < 1.0e-8
	@test norm(ContactControl.B_fast(model, q0s) - ContactControl.B_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactControl.A_fast(model, q0s) - ContactControl.A_func(model, q0s), Inf) < 1.0e-8
    @test norm(ContactControl.C_fast(model, q0s, q̇0s) - ContactControl.C_func(model, q0s, q̇0s), Inf) < 1.0e-8
end

@testset "Model Methods: Fast Dynamics (Particle)" begin
	s = ContactControl.get_simulation("particle", "flat_3D_lc", "flat")
    model = s.model
	env = s.env

	# Setup variables
    T = Float64
    nq = model.dim.q
    nu = model.dim.u
    nw = model.dim.w
    nc = model.dim.c
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
    @test norm(ContactControl.d_fast(model, hs, q0s, q1s, u1s, w1s, λ1s, q2s) -
        ContactControl.dynamics(model, hs, q0s, q1s, u1s, w1s, λ1s, q2s), Inf) < 1.0e-8
end
