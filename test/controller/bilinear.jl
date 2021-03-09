################################################################################
# Test
################################################################################

@testset "Bilinear Approximations" begin
	names = ["particle", "quadruped"]
	for name in names
		model = ContactControl.get_model(name)
		nz = ContactControl.num_var(model)
		nθ = ContactControl.num_data(model)
		z = rand(SizedVector{nz,Float64})
		θ = rand(SizedVector{nθ,Float64})
		κ = 1e-5
		r0 = rand(SizedVector{nz,Float64})
		r1 = rand(SizedVector{nz,Float64})
		rz0 = spzeros(nz,nz)
		rz0 = similar(model.spa.rz_sp, Float64)
		rz1 = deepcopy(rz0)
		rθ = zeros(nz,nθ)
		lin = ContactControl.LinStep14(model, z, θ, κ)

		# Test r!
		model.res.r(r0, z, θ, κ)
		ContactControl.r_approx!(model, lin, r1, z, θ, κ)
		@test norm(r0 - r1, Inf) < 1e-8

		α = 1.1
		ContactControl.r_approx!(model, lin, r1, α*z, α*θ, κ)
		@test r0 != r1

		# Test rz!
		function rz_approx_FD!(model::ContactDynamicsModel, lin::LinStep14, rz::AbstractMatrix{T},
			z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
			nz = ContactControl.num_var(model)
			function f(z)
				r = zeros(SizedVector{nz,eltype(z)})
				ContactControl.r_approx!(model, lin, r, z, θ, κ)
				return Vector(r)
			end
			return rz = ForwardDiff.jacobian(f, Vector(z))
		end

		model.res.rz(rz0, z, θ, κ)
		ContactControl.rz_approx!(model, lin, rz1, z, θ, κ)
		rz_FD = rz_approx_FD!(model, lin, rz1, z, θ, κ)
		@test norm(rz0 - rz1, Inf) < 1e-8
		@test norm(rz1 - rz_FD, Inf) < 1e-8

		α = 1.1
		ContactControl.rz_approx!(model, lin, rz1, α*z, α*θ, κ)
		rz_FD = rz_approx_FD!(model, lin, rz1, α*z, α*θ, κ)
		@test rz0 != rz1
		@test norm(rz_FD - rz1, Inf) < 1e-8
	end
end
