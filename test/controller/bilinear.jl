################################################################################
# Test
################################################################################

@testset "Bilinear Approximations" begin
	names = ["particle", "quadruped"]
	for name in names
		model = get_model(name)
		nz = num_var(model)
		nθ = num_data(model)
		z_ = rand(SizedVector{nz,Float64})
		θ_ = rand(SizedVector{nθ,Float64})
		κ_ = 1e-5
		r0_ = rand(SizedVector{nz,Float64})
		r1_ = rand(SizedVector{nz,Float64})
		rz0_ = spzeros(nz,nz)
		rz0_ = similar(model.spa.rz_sp, Float64)
		rz1_ = deepcopy(rz0_)
		rθ_ = zeros(nz,nθ)
		lin = LinStep14(model, z_, θ_, κ_)

		# Test r!
		model.res.r(r0_, z_, θ_, κ_)
		r_approx!(model, lin, r1_, z_, θ_, κ_)
		@test norm(r0_ - r1_, Inf) < 1e-8

		α = 1.1
		r_approx!(model, lin, r1_, α*z_, α*θ_, κ_)
		@test r0_ != r1_

		# Test rz!
		function rz_approx_FD!(model::ContactDynamicsModel, lin::LinStep14, rz::AbstractMatrix{T},
			z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
			nz = num_var(model)
			function f(z)
				r = zeros(SizedVector{nz,eltype(z)})
				r_approx!(model, lin, r, z, θ, κ)
				return r
			end
			return rz = ForwardDiff.jacobian(f, z)
		end

		model.res.rz(rz0_, z_, θ_, κ_)
		rz_approx!(model, lin, rz1_, z_, θ_, κ_)
		rz_FD = rz_approx_FD!(model, lin, rz1_, z_, θ_, κ_)
		@test norm(rz0_ - rz1_, Inf) < 1e-8
		@test norm(rz1_ - rz_FD, Inf) < 1e-8

		α = 1.1
		rz_approx!(model, lin, rz1_, α*z_, α*θ_, κ_)
		rz_FD = rz_approx_FD!(model, lin, rz1_, α*z_, α*θ_, κ_)
		@test rz0_ != rz1_
		@test norm(rz_FD - rz1_, Inf) < 1e-8
	end
end
