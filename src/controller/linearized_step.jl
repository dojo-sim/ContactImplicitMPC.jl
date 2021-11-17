# struct LinearizedStep{T}
# 	z::Vector{T}
# 	θ::Vector{T}
# 	κ::Vector{T}
# 	r::Vector{T}
# 	rz::Matrix{T}
# 	rθ::Matrix{T}
# end

# function LinearizedStep(s::Simulation, z::Vector{T}, θ::Vector{T}, κ::T) where T
# 	model = s.model
# 	env = s.env

# 	nz = num_var(model, env)
# 	nθ = num_data(model)

# 	z0 = zero(z)
# 	θ0 = zero(θ)
# 	κ0 = [κ]
# 	r0 = zeros(nz)
# 	rz0 = zeros(nz, nz)
# 	rθ0 = zeros(nz, nθ)

# 	s.res.r!(r0, z0, θ0, κ0)
# 	s.res.rz!(rz0, z0, θ0)
# 	s.res.rθ!(rθ0, z0, θ0)

# 	return LinearizedStep(z0, θ0, κ0, r0, rz0, rθ0)
# end

# """
# 	Create LinearizedStep.
# """
# function LinearizedStep(s::Simulation)
# 	model = s.model
# 	env = s.env

# 	nz = num_var(model, env)
# 	nθ = num_data(model)
# 	z0 = zeros(nz)
# 	θ0 = zeros(nθ)
# 	κ0 = [0.0]

# 	return LinearizedStep(s, z0, θ0, κ0)
# end

# function update!(lin::LinearizedStep, s::Simulation, z, θ)
# 	lin.z .= z
# 	lin.θ .= θ
# 	s.res.r!(lin.r, z, θ, lin.κ)
# 	s.res.rz!(lin.rz, z, θ)
# 	s.res.rθ!(lin.rθ, z, θ)
# 	return nothing
# end

mutable struct LinearizedStep{T}
	z::AbstractVector{T}
	θ::AbstractVector{T}
	κ::T
	r::AbstractVector{T}
	rz::AbstractMatrix{T}
	rθ::AbstractMatrix{T}
end

function LinearizedStep(s::Simulation, z::AbstractVector{T}, θ::AbstractVector{T}, κ::T) where T
	model = s.model
	env = s.env

	nz = num_var(model, env)
	nθ = num_data(model)

	z0 = SizedVector{nz,T}(z)
	θ0 = SizedVector{nθ,T}(θ)
	κ0 = κ
	r0 = zeros(SizedVector{nz,T})
	rz0 = zeros(nz, nz)
	rθ0 = zeros(nz, nθ)

	s.res.r!(r0, z0, θ0, κ0)
	s.res.rz!(rz0, z0, θ0)
	s.res.rθ!(rθ0, z0, θ0)

	return LinearizedStep(z0, θ0, κ0, r0, rz0, rθ0)
end

"""
	Create LinearizedStep.
"""
function LinearizedStep(s::Simulation)
	model = s.model
	env = s.env

	nz = num_var(model, env)
	nθ = num_data(model)
	z0 = zeros(SizedVector{nz})
	θ0 = zeros(SizedVector{nθ})
	κ0 = 0.0

	return LinearizedStep(s, z0, θ0, κ0)
end

function update!(lin::LinearizedStep, s::Simulation, z, θ)
	lin.z .= z
	lin.θ .= θ
	s.res.r!(lin.r, z, θ, lin.κ)
	s.res.rz!(lin.rz, z, θ)
	s.res.rθ!(lin.rθ, z, θ)
	return nothing
end