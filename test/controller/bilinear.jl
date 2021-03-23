################################################################################
# Test
################################################################################

@testset "Controller: Bilinear Approximations" begin
	names = ["particle", "quadruped"]

	for name in names
		model = ContactControl.get_model(name)

		nz = ContactControl.num_var(model)
		nθ = ContactControl.num_data(model)
		z = rand(nz)
		θ = rand(nθ)
		κ = 1e-5

		r0 = rand(nz)
		r1 = rand(nz)
		r2 = rand(nz)
		rz0 = spzeros(nz,nz)
		rz0 = similar(model.spa.rz_sp, Float64)
		rz1 = deepcopy(rz0)
		rθ0 = zeros(nz,nθ)
		lin = ContactControl.LinStep(model, z, θ, κ)
		@test norm(lin.z0 - z) < 1e-8
		@test norm(lin.θ0 - θ) < 1e-8
		@test norm(lin.κ0 - κ) < 1e-8

		# Test r!: the linear approximaton is exact at the linearization point
		model.res.r(r0, z, θ, κ)
		ContactControl.r_approx!(lin, r1, z, θ, κ)
		lin.methods.r!(r2, z, θ, κ)
		@test norm(r0 - r1, Inf) < 1e-8
		@test norm(r0 - r2, Inf) < 1e-8

		# Test rz!: the linear approximaton is exact at the linearization point
		function rz_approx_FD!(model::ContactDynamicsModel, lin::LinStep, rz::AbstractMatrix{T},
			z::AbstractVector{T}, θ::AbstractVector{T}) where {T}
			nz = ContactControl.num_var(model)
			function f(z)
				r = zeros(SizedVector{nz,eltype(z)})
				ContactControl.r_approx!(lin, r, z, θ, κ)
				return Vector(r)
			end
			return rz = ForwardDiff.jacobian(f, Vector(z))
		end
		model.res.rz(rz0, z, θ)
		ContactControl.rz_approx!(lin, rz1, z, θ)
		rz_FD = rz_approx_FD!(model, lin, rz1, z, θ)
		@test norm(rz0 - rz1, Inf) < 1e-8
		@test norm(rz1 - rz_FD, Inf) < 1e-8

		# Test r!: at a point away from the linearization point
		α = 2.0
		β = 1.1
		r1 = rand(nz)
		r2 = rand(nz)
		r3 = rand(nz)
		r4 = rand(nz)
		model.res.r(r1, α*z, β*θ, κ)
		ContactControl.r_approx!(lin, r2, α*z, β*θ, κ)
		model.approx.r(r3, α*z, β*θ, κ, lin.z0, lin.θ0, lin.r0, lin.rz0, lin.rθ0)
		lin.methods.r!(r4, α*z, β*θ, κ)
		@test norm(r2 - r3, Inf) < 1e-8
		@test norm(r2 - r4, Inf) < 1e-8
		norm(r1 - r3, Inf)

		# Test rz!: at a point away from the linearization point
		rz1 = rand(nz,nz)
		rz2 = rand(nz,nz)
		rz3 = rand(nz,nz)
		rz4 = rand(nz,nz)
		model.res.rz(rz1, α*z, β*θ)
		ContactControl.rz_approx!(lin, rz2, α*z, β*θ)
		model.approx.rz(rz3, α*z, lin.rz0)
		lin.methods.rz!(rz4, α*z, β*θ)
		@test norm(rz2 - rz3, Inf) < 1e-8
		@test norm(rz2 - rz4, Inf) < 1e-8

		# Test rθ!: at a point away from the linearization point
		rθ1 = rand(nz,nθ)
		rθ2 = rand(nz,nθ)
		rθ3 = rand(nz,nθ)
		rθ4 = rand(nz,nθ)
		zeros(nz,nθ)
		model.res.rθ(rθ1, α*z, β*θ)
		ContactControl.rθ_approx!(lin, rθ2, α*z, β*θ)
		model.approx.rθ(rθ3, lin.rθ0)
		lin.methods.rθ!(rθ4, α*z, β*θ)
		@test norm(rθ2 - rθ3, Inf) < 1e-8
		@test norm(rθ2 - rθ4, Inf) < 1e-8
	end
end
