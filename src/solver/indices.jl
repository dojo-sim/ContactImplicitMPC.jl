################################################################################
# Indices
################################################################################

mutable struct IndicesOptimization 
	# Set the residual to 0
	# r(z) = 0
	# z <- z + Δz

	# Dimensions
	nz::Int # dimension of the optimization variable z
	nΔ::Int # dimension of the optimization variable Δz and of the residual r
	ny::Int # dimension of the bilinear vaiables
	nquat::Int # number of quaternions

	# Variables
	dynz::Vector{Int} # indices of the variables associated with the dynamics constraints in z
	dynΔ::Vector{Int} # indices of the variables associated with the dynamics constraints in Δz
	ortz::Vector{Vector{Int}} # indices of the variables associated with the positive ORThant constraints in z
	ortΔ::Vector{Vector{Int}} # indices of the variables associated with the positive ORThant constraints in Δz
	socz::Vector{Vector{Vector{Int}}} # indices of the variables associated with the Second Order Cone constraints in z
	socΔ::Vector{Vector{Vector{Int}}} # indices of the variables associated with the Second Order Cone constraints in Δz

	# Residual
	equr::Vector{Int} # indices of the residual associated with the EQUality constraints in r
	ortr::Vector{Int} # indices of the residual associated with the positive ORThant constraints in r
	socr::Vector{Int} # indices of the residual associated with the Second Order Cone constraints in r
	dyn::Vector{Int} # indices of the residual associated with the dynamics constraints in r
	rst::Vector{Int} # indices of the residual associated with the remaining constraints in r
	bil::Vector{Int} # indices of the residual associated with the bilinear constraints in r
	alt::Vector{Int} # indices of the residual associated with the altitude constraints in r 
end

function IndicesOptimization()
	v1 = Vector{Int}()
	v2 = Vector{Vector{Int}}()
	v3 = Vector{Vector{Vector{Int}}}()

	s = IndicesOptimization(
		0, 0, 0, 0,
		v1, v1, v2, v2, v3, v3,
		v1, v1, v2, v1, v1, v1, v1)
	return s
end
