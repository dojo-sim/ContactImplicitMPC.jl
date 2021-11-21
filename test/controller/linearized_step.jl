@testset "Controller: Linearized Approximations" begin
	names = ["particle", "quadruped"]

	ss = [get_simulation("particle", "quadratic_bowl_3D_lc", "quadratic"),
	      get_simulation("quadruped", "flat_2D_lc", "flat")]

	for s in ss
		nz = ContactImplicitMPC.num_var(s.model, s.env)
		nθ = ContactImplicitMPC.num_data(s.model)
		z = rand(nz)
		θ = rand(nθ)
		κ = 1.0e-5

		r0 = rand(nz)
		rz0 = similar(s.rz, Float64)

		lin = ContactImplicitMPC.LinearizedStep(s, z, θ, κ)

		@test norm(lin.z - z) < 1.0e-8
		@test norm(lin.θ - θ) < 1.0e-8
		@test norm(lin.κ[1] - κ) < 1.0e-8

		# Test r!: the linear linearization is exact at the linearization point
		s.res.r!(r0, z, θ, [κ])
		@test norm(r0 - lin.r, Inf) < 1.0e-8

		s.res.rz!(rz0, z, θ)
		@test norm(rz0 - lin.rz, Inf) < 1.0e-8

		#TODO: add more tests
	end
end
