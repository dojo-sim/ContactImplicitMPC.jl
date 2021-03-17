mutable struct LinStep{T}
	z0::AbstractVector{T}
	θ0::AbstractVector{T}
	κ0::T
	r0::AbstractVector{T}
	rz0::AbstractMatrix{T}
	rθ0::AbstractMatrix{T}
	bil_terms::Any
	bil_vars::Any
end

function LinStep(model::ContactDynamicsModel, z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
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
	model.res.rz(rz0, z0, θ0)
	model.res.rθ(rθ0, z0, θ0)
	bil_terms, bil_vars = get_bilinear_indices(model)
	return LinStep{T}(z0, θ0, κ0, r0, rz0, rθ0, bil_terms, bil_vars)
end

"""
	Create dummy LinStep.
"""
function LinStep(model::ContactDynamicsModel)
	nz = num_var(model)
	nθ = num_data(model)
	z0 = zeros(SizedVector{nz})
	θ0 = zeros(SizedVector{nθ})
	κ0 = 0.0
	return LinStep(model, z0, θ0, κ0)
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

function bil_addition!(out::AbstractVector, i::SVector{n,Int}, a::AbstractVector,
	b::AbstractVector, ρ::T) where {n,T}
	out[i] = a.*b .- ρ
	return nothing
end

"""
	r_approx!(lin::LinStep, r::AbstractVector{T1},
	z::AbstractVector{T1}, θ::AbstractVector{T}, κ::T) where {T,T1}
Compute an approximate residual. The approximation results from the linearization of the non-linear
terms in the residual about a reference point. The bilinear terms (complementarity constraints) are
not linearized.
"""
function r_approx!(lin::LinStep, r::AbstractVector{T1},
	z::AbstractVector{T1}, θ::AbstractVector{T}, κ::T) where {T,T1}
	@assert norm(κ - lin.κ0)/κ < 1e-10
	r .= lin.r0 + lin.rz0 * (z-lin.z0) + lin.rθ0 * (θ-lin.θ0)
	for i = 1:length(lin.bil_terms)
		t = lin.bil_terms[i]
		v1 = lin.bil_vars[i][1]
		v2 = lin.bil_vars[i][2]
		# r[t] = z[v1].*z[v2] .- κ
		bil_addition!(r, t, z[v1], z[v2], κ)
	end
    return nothing
end

"""
	rz_approx!(lin::LinStep, rz::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
Compute an approximate residual jacobian with respect to z. The approximation results from the linearization of the non-linear
terms in the residual about a reference point. The bilinear terms (complementarity constraints) are
not linearized.
"""
function rz_approx!(lin::LinStep, rz::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}) where {T}
	rz .= lin.rz0
	for i = 1:length(lin.bil_terms)
		t = lin.bil_terms[i]
		v1 = lin.bil_vars[i][1]
		v2 = lin.bil_vars[i][2]
		rz[t,v1] .= Diagonal(z[v2])
		rz[t,v2] .= Diagonal(z[v1])
	end
    return nothing
end

"""
	rθ_approx!(lin::LinStep, rz::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
Compute an approximate residual jacobian with respect to θ The approximation results from the linearization of the non-linear
terms in the residual about a reference point. The bilinear terms (complementarity constraints) are
not linearized.
"""
function rθ_approx!(lin::LinStep, rθ::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}) where {T}
	rθ .= lin.rθ0
    return nothing
end
