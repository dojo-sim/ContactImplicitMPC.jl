mutable struct LinStep14{T}
	z0::AbstractVector{T}
	θ0::AbstractVector{T}
	κ0::T
	r0::AbstractVector{T}
	rz0::AbstractMatrix{T}
	rθ0::AbstractMatrix{T}
	bil_terms::Any
	bil_vars::Any
end

function LinStep14(model::ContactDynamicsModel, z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
	nz = num_var(model)
	nθ = num_data(model)
	z0 = SizedVector{nz,T}(z)
	θ0 = SizedVector{nθ,T}(θ)
	κ0 = κ
	r0 = zeros(SizedVector{nz,T})
	rz0 = spzeros(nz,nz)
	rz0 = similar(model.spa.rz_sp, T)
	rθ0 = zeros(nz, nθ)
	model.res.r(r0, z0, θ0, κ0)
	model.res.rz(rz0, z0, θ0, κ0)
	model.res.rθ(rθ0, z0, θ0, κ0)
	bil_terms, bil_vars = get_bilinear_indices(model)
	return LinStep14{T}(z0, θ0, κ0, r0, rz0, rθ0, bil_terms, bil_vars)
end

function get_bilinear_indices(model::ContactDynamicsModel)
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b

	# lin_terms = [Vector(1:nq),
	# 			 Vector(nq .+ (1:nc)),
	# 			 Vector(nq+2nc .+ (1:nb)),
	# 			 Vector(nq+2nc+nb .+ (1:nc))]
	bil_terms = [SVector{nc,Int}(nq+nc .+ (1:nc)),
				 SVector{nc,Int}(nq+3nc+nb .+ (1:nc)),
				 SVector{nb,Int}(nq+4nc+nb .+ (1:nb))]
	bil_vars = [[SVector{nc,Int}(nq .+ (1:nc)), SVector{nc}(nq+2nc+2nb .+ (1:nc))],  # γ1, s1
				[SVector{nc,Int}(nq+nc+nb .+ (1:nc)), SVector{nc}(nq+3nc+2nb .+ (1:nc))],  # ψ, s2
				[SVector{nb,Int}(nq+nc .+ (1:nb)), SVector{nb}(nq+2nc+nb .+ (1:nb))]] # b2, η
	return bil_terms, bil_vars
end

function bil_addition!(out::AbstractVector{T}, i::SVector{n,Int}, a::SizedVector{n,T},
	b::SizedVector{n,T}, ρ::T) where {n,T}
	out[i] = a.*b .- ρ
	return nothing
end
