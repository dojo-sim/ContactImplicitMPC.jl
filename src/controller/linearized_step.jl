mutable struct LinearizedStep{T}
	z::AbstractVector{T}
	θ::AbstractVector{T}
	κ::T
	r::AbstractVector{T}
	rz::AbstractMatrix{T}
	rθ::AbstractMatrix{T}
	terms::Any
	vars::Any
	methods::ResidualMethods
end

function LinearizedStep(model::ContactDynamicsModel, z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where T
	nz = num_var(model)
	nθ = num_data(model)

	z0 = SizedVector{nz,T}(z)
	θ0 = SizedVector{nθ,T}(θ)
	κ0 = κ
	r0 = zeros(SizedVector{nz,T})
	rz0 = zeros(nz, nz)
	rθ0 = zeros(nz, nθ)

	model.res.r!(r0, z0, θ0, κ0)
	model.res.rz!(rz0, z0, θ0)
	model.res.rθ!(rθ0, z0, θ0)

	terms, vars = get_bilinear_indices(model)

	function r_linearized!(r, z, θ, κ)
		model.linearized.r!(r, z, θ, κ, z0, θ0, r0, rz0, rθ0)
		return nothing
	end

	function rz_linearized!(rz, z, θ)
		model.linearized.rz!(rz, z, rz0)
		return nothing
	end

	function rθ_linearized!(rθ, z, θ)
		model.linearized.rθ!(rθ, rθ0)
		return nothing
	end

	methods = ResidualMethods(r_linearized!, rz_linearized!, rθ_linearized!)

	return LinearizedStep(z0, θ0, κ0, r0, rz0, rθ0, terms, vars, methods)
end

"""
	Create LinearizedStep.
"""
function LinearizedStep(model::ContactDynamicsModel)
	nz = num_var(model)
	nθ = num_data(model)
	z0 = zeros(SizedVector{nz})
	θ0 = zeros(SizedVector{nθ})
	κ0 = 0.0
	return LinearizedStep(model, z0, θ0, κ0)
end

function get_bilinear_indices(model::ContactDynamicsModel)
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b

	terms = [SVector{nc,Int}(nq + nc .+ (1:nc)),
				 SVector{nc,Int}(nq + 3nc + nb .+ (1:nc)),
				 SVector{nb,Int}(nq + 4nc + nb .+ (1:nb))]

	vars = [[SVector{nc,Int}(nq .+ (1:nc)), SVector{nc}(nq + 2nc + 2nb .+ (1:nc))],  # γ1, s1
				[SVector{nc,Int}(nq + nc + nb .+ (1:nc)), SVector{nc}(nq + 3nc + 2nb .+ (1:nc))],  # ψ, s2
				[SVector{nb,Int}(nq + nc .+ (1:nb)), SVector{nb}(nq + 2nc + nb .+ (1:nb))]] # b1, η

	return terms, vars
end

"""
	r_linearized!(lin::LinearizedStep, r::AbstractVector{T1},
	z::AbstractVector{T1}, θ::AbstractVector{T}, κ::T) where {T,T1}
Compute an linearizedimate residual. The linearizedimation results from the linearization of the non-linear
terms in the residual about a reference point. The bilinear terms (complementarity constraints) are
not linearized.
"""
function r_linearized!(lin::LinearizedStep, r::AbstractVector{T1},
	z::AbstractVector{T1}, θ::AbstractVector{T}, κ::T) where {T,T1}

	@assert norm(κ - lin.κ) / κ < 1e-10

	r .= lin.r + lin.rz * (z - lin.z) + lin.rθ * (θ - lin.θ)

	for i = 1:length(lin.terms)
		t = lin.terms[i]
		v1 = lin.vars[i][1]
		v2 = lin.vars[i][2]
		r[t] = z[v1] .* z[v2] .- κ
	end

    return nothing
end

"""
	rz_linearized!(lin::LinearizedStep, rz::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
Compute an linearizedimate residual jacobian with respect to z. The linearizedimation results from the linearization of the non-linear
terms in the residual about a reference point. The bilinear terms (complementarity constraints) are
not linearized.
"""
function rz_linearized!(lin::LinearizedStep, rz::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}) where {T}

	rz .= lin.rz

	for i = 1:length(lin.terms)
		t = lin.terms[i]
		v1 = lin.vars[i][1]
		v2 = lin.vars[i][2]
		rz[t,v1] .= Diagonal(z[v2])
		rz[t,v2] .= Diagonal(z[v1])
	end

    return nothing
end

"""
	rθ_linearized!(lin::LinearizedStep, rz::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where {T}
Compute an linearizedimate residual jacobian with respect to θ The linearizedimation results from the linearization of the non-linear
terms in the residual about a reference point. The bilinear terms (complementarity constraints) are
not linearized.
"""
function rθ_linearized!(lin::LinearizedStep, rθ::AbstractMatrix{T},
	z::AbstractVector{T}, θ::AbstractVector{T}) where {T}

	rθ .= lin.rθ

    return nothing
end
