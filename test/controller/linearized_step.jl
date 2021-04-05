################################################################################
# Test
################################################################################

@testset "Controller: Linearized Approximations" begin
	names = ["particle", "quadruped"]

	for name in names
		model = ContactControl.get_model(name)

		nz = ContactControl.num_var(model)
		nθ = ContactControl.num_data(model)
		z = rand(nz)
		θ = rand(nθ)
		κ = 1.0e-5

		r0 = rand(nz)
		rz0 = similar(model.spa.rz_sp, Float64)

		lin = ContactControl.LinearizedStep(model, z, θ, κ)

		@test norm(lin.z - z) < 1.0e-8
		@test norm(lin.θ - θ) < 1.0e-8
		@test norm(lin.κ - κ) < 1.0e-8

		# Test r!: the linear linearization is exact at the linearization point
		model.res.r!(r0, z, θ, κ, nothing)
		@test norm(r0 - lin.r, Inf) < 1.0e-8

		model.res.rz!(rz0, z, θ, nothing)
		@test norm(rz0 - lin.rz, Inf) < 1.0e-8

		#TODO: add more tests
	end
end
