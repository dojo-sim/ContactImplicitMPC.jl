mutable struct LinearizedStep{T}
	z::AbstractVector{T}
	θ::AbstractVector{T}
	κ::T
	r::AbstractVector{T}
	rz::AbstractMatrix{T}
	rθ::AbstractMatrix{T}
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
	return LinearizedStep(z0, θ0, κ0, r0, rz0, rθ0)
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

	terms = [SVector{nc,Int}(nq + nb + 2nc .+ (1:nc)),
			 SVector{nb,Int}(nq + nb + 3nc .+ (1:nb)),
			 SVector{nc,Int}(nq + nb + 3nc + nb .+ (1:nc))]

	vars = [[SVector{nc,Int}(nq .+ (1:nc)),           SVector{nc}(nq + 2nc + 2nb .+ (1:nc))],  # γ1, s1
			[SVector{nb,Int}(nq + nc .+ (1:nb)),      SVector{nb}(nq + 2nc +  nb .+ (1:nb))], # b1, η
			[SVector{nc,Int}(nq + nc + nb .+ (1:nc)), SVector{nc}(nq + 3nc + 2nb .+ (1:nc))],  # ψ, s2
			]
	return terms, vars
end

function update!(lin::LinearizedStep, model::ContactDynamicsModel, z, θ)
	lin.z .= z
	lin.θ .= θ
	model.res.r!(lin.r, z, θ, lin.κ)
	model.res.rz!(lin.rz, z, θ)
	model.res.rθ!(lin.rθ, z, θ)
	return nothing
end
