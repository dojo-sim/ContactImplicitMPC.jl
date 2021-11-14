# TODO
# we need to have nquat directly inside the model

################################################################################
# z indices
################################################################################

"""
	Returns the indices of q2 in z and Δq2 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_q2(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nquat = 0
	iq2 = Vector(1:nq - nquat)
	return iq2
end

"""
	Returns the indices of γ1 in z and Δγ1 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_γ1(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	off = (nq - nquat)
	iγ1 = Vector(off .+ (1:nc))
	return iγ1
end

"""
	Returns the indices of b1 in z and Δb1 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_b1(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc
	ib1 = Vector(off .+ (1:nb))
	return ib1
end

"""
	Returns the indices of ψ1 in z and Δψ1 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_ψ1(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb
	iψ1 = Vector(off .+ (1:nc))
	return iψ1
end

"""
	Returns the indices of s1 in z and Δs1 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_s1(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc
	is1 = Vector(off .+ (1:nc))
	return is1
end

"""
	Returns the indices of η1 in z and Δη1 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_η1(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc + nc
	iη1 = Vector(off .+ (1:nb))
	return iη1
end

"""
	Returns the indices of s2 in z and Δs2 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_s2(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc + nc + nb
	is2 = Vector(off .+ (1:nc))
	return is2
end


################################################################################
# θ indices
################################################################################

"""
	Returns the indices of q0 in θ.
"""
function index_q0(model::Model)
	nq = model.nq
	iq0 = Vector(1:nq)
	return iq0
end

"""
	Returns the indices of q1 in θ.
"""
function index_q1(model::Model)
	nq = model.nq
	off = nq
	iq1 = Vector(off .+ (1:nq))
	return iq1
end

"""
	Returns the indices of u1 in θ.
"""
function index_u1(model::Model)
	nq = model.nq
	nu = model.nu
	off = nq + nq
	iu1 = Vector(off .+ (1:nu))
	return iu1
end

"""
	Returns the indices of w1 in θ.
"""
function index_w1(model::Model)
	nq = model.nq
	nu = model.nu
	nw = model.nw
	off = nq + nq + nu
	iw1 = Vector(off .+ (1:nw))
	return iw1
end

"""
	Returns the indices of μ in θ.
"""
function index_μ(model::Model)
	nq = model.nq
	nu = model.nu
	nw = model.nw
	off = nq + nq + nu + nw
	iμ = Vector(off .+ (1:1))
	return iμ
end

"""
	Returns the indices of h in θ.
"""
function index_h(model::Model)
	nq = model.nq
	nu = model.nu
	nw = model.nw
	off = nq + nq + nu + nw + 1
	ih = Vector(off .+ (1:1))
	return ih
end

################################################################################
# Residual indices
################################################################################

"""
	Returns the indices of the dynamics equation in the residual r.
"""
function index_dyn(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nquat = 0
	idyn = Vector(1:nq - nquat)
	return idyn
end

"""
	Returns the indices of the impact equation in the residual r.
"""
function index_imp(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	off = (nq - nquat)
	iimp = Vector(off .+ (1:nc))
	return iimp
end

"""
	Returns the indices of the maximum dissipation principle equation in the residual r.
"""
function index_mdp(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc
	imdp = Vector(off .+ (1:nb))
	return imdp
end

"""
	Returns the indices of the friction cone equation in the residual r.
"""
function index_fri(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb
	ifri = Vector(off .+ (1:nc))
	return ifri
end

"""
	Returns the indices of the bilinear impact equation in the residual r.
"""
function index_bimp(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc
	ibimp = Vector(off .+ (1:nc))
	return ibimp
end

"""
	Returns the indices of the bilinear maximum dissipation principle equation in the residual r.
"""
function index_bmdp(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc + nc
	ibmdp = Vector(off .+ (1:nb))
	return ibmdp
end

"""
	Returns the indices of the bilinear friction equation in the residual r.
"""
function index_bfri(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc + nc + nb
	ibfri = Vector(off .+ (1:nc))
	return ibfri
end

################################################################################
# Aggregated Indices
################################################################################

"""
	Returns the indices of z and Δz, where 3 groups are formed.

	w1 = (q2)
	w2 = (γ1, b1, ψ1)
	w3 = (s1, η1, s2)

	Δw1 = (Δq2)
	Δw2 = (Δγ1, Δb1, Δψ1)
	Δw3 = (Δs1, Δη1, Δs2)

	Note that Δq2 does not have the same size as q2 when optimizing over
	non-euclidean spaces.
"""
function linearization_var_index(model::Model, env::Environment; quat::Bool = false)
	iq2 = index_q2(model, env, quat = quat)
	iγ1 = index_γ1(model, env, quat = quat)
	ib1 = index_b1(model, env, quat = quat)
	iψ1 = index_ψ1(model, env, quat = quat)
	is1 = index_s1(model, env, quat = quat)
	iη1 = index_η1(model, env, quat = quat)
	is2 = index_s2(model, env, quat = quat)
	iw1 = iq2
	iw2 = [iγ1; ib1; iψ1]
	iw3 = [is1; iη1; is2]
	return iw1, iw2, iw3
end

"""
	Returns the indices of the residual r, where 3 groups are formed.
	r = (rdyn, rrst, rbil)
"""
function linearization_term_index(model::Model, env::Environment; quat::Bool = false)
	# dyn = [dyn]
	# rst = [s1  - ..., ≡ ialt
	#        η1  - ...,
	#        s2  - ...,]
	# bil = [γ1 .* s1 .- κ;
	#        b1 .* η1 .- κ;
	#        ψ1 .* s2 .- κ]
	idyn = index_dyn(model, env, quat = quat)
	iimp = index_imp(model, env, quat = quat)
	imdp = index_mdp(model, env, quat = quat)
	ifri = index_fri(model, env, quat = quat)
	ibimp = index_bimp(model, env, quat = quat)
	ibmdp = index_bmdp(model, env, quat = quat)
	ibfri = index_bfri(model, env, quat = quat)

	irst = [iimp; imdp; ifri]
	ibil = [ibimp; ibmdp; ibfri]
	ialt = iimp
	return idyn, irst, ibil, ialt
end

"""
	Returns the positive orthant indices in z or Δz.
"""
function index_ort(model::Model, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
	iγ1 = index_γ1(model, env, quat = quat)
	ib1 = index_b1(model, env, quat = quat)
	iψ1 = index_ψ1(model, env, quat = quat)
	is1 = index_s1(model, env, quat = quat)
	iη1 = index_η1(model, env, quat = quat)
	is2 = index_s2(model, env, quat = quat)
	return [[iγ1; ib1; iψ1], [is1; iη1; is2]]
end

"""
	Returns the positive orthant indices in z or Δz.
"""
function index_ort(model::Model, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	iγ1 = index_γ1(model, env, quat = quat)
	is1 = index_s1(model, env, quat = quat)
	return [iγ1, is1]
end

"""
	Returns the second order cone indices in z or Δz.
"""
index_soc(model::Model, env::Environment{<:World,LinearizedCone}; quat::Bool = false) = Vector{Vector{Int}}[]

"""
	Returns the second order cone indices in z or Δz.
"""
function index_soc(model::Model, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	nc = model.nc
	nf = friction_dim(env)

	ib1 = index_b1(model, env, quat = quat) # primal cones: vector part
	iψ1 = index_ψ1(model, env, quat = quat) # primal cones: scalar part
	iη1 = index_η1(model, env, quat = quat) # dual cones: vector part
	is2 = index_s2(model, env, quat = quat) # dual cones: scalar part

	[[[iψ1[i]; iη1[(i - 1) * nf .+ (1:nf)]], [is2[i]; ib1[(i - 1) * nf .+ (1:nf)]]] for i = 1:nc]
end

function num_var(model::Model, env::Environment; quat::Bool = false)
	nq = model.nq
	nc = model.nc
	nquat = 0
	nb = nc * friction_dim(env)
	(nq - nquat) + nc + nb + nc + nc + nb + nc
end

function num_data(model::Model)
	nq = model.nq
	nu = model.nu
	nw = model.nw
	nq + nq + nu + nw + 1 + 1
end

function num_bilinear(model::Model, env::Environment)
	nc = model.nc
	nb = nc * friction_dim(env)
	nc + nb + nc
end

################################################################################
# Packing and Unpacking Variables
################################################################################

function unpack_θ(model::Model, θ)
	iq0 = index_q0(model)
	iq1 = index_q1(model)
	iu1 = index_u1(model)
	iw1 = index_w1(model)
	iμ  = index_μ(model)
	ih  = index_h(model)

	q0 = θ[iq0]
	q1 = θ[iq1]
	u1 = θ[iu1]
	w1 = θ[iw1]
	μ = θ[iμ]
	h = θ[ih]
	return q0, q1, u1, w1, μ, h
end

function pack_θ(model::Model, q0, q1, u1, w1, μ, h)
	return [q0; q1; u1; w1; μ; h]
end

function unpack_z(model::Model, env::Environment, z)
	iq2 = index_q2(model, env, quat = false)
	iγ1 = index_γ1(model, env, quat = false)
	ib1 = index_b1(model, env, quat = false)
	iψ1 = index_ψ1(model, env, quat = false)
	is1 = index_s1(model, env, quat = false)
	iη1 = index_η1(model, env, quat = false)
	is2 = index_s2(model, env, quat = false)

	# system variables
	q2 = z[iq2]
	γ1 = z[iγ1]
	b1 = z[ib1]
	ψ1 = z[iψ1]
	s1 = z[is1]
	η1 = z[iη1]
	s2 = z[is2]
	return q2, γ1, b1, ψ1, s1, η1, s2
end

function pack_z(model::Model, env::Environment{<:World,LinearizedCone}, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_func(model, env, q2)
	s2 = model.μ_world * γ1 .- E_func(model, env) * b1
	return pack_z(model, env, q2, γ1, b1, ψ1, s1, η1, s2)
end

function pack_z(model::Model, env::Environment{<:World,NonlinearCone}, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_func(model, env, q2)
	s2 = model.μ_world .* γ1
	return pack_z(model, env, q2, γ1, b1, ψ1, s1, η1, s2)
end

function pack_z(model::Model, env::Environment, q2, γ1, b1, ψ1, s1, η1, s2)
	return [q2; γ1; b1; ψ1; s1; η1; s2]
end

################################################################################
# Aggregated indices
################################################################################

function index_equr(model::Model, env::Environment; quat::Bool = false)
	idyn = index_dyn(model, env, quat = quat)
	iimp = index_imp(model, env, quat = quat)
	imdp = index_mdp(model, env, quat = quat)
	ifri = index_fri(model, env, quat = quat)
	iequ = [idyn; iimp; imdp; ifri]
	return iequ
end

function index_ortr(model::Model, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
	ibimp = index_bimp(model, env, quat = quat)
	ibmdp = index_bmdp(model, env, quat = quat)
	ibfri = index_bfri(model, env, quat = quat)
	iortr = [ibimp; ibmdp; ibfri]
	return iortr
end

function index_ortr(model::Model, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	ibimp = index_bimp(model, env, quat = quat)
	iortr = ibimp
	return iortr
end

function index_socr(model::Model, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
	isocr = Vector{Int64}()
	return isocr
end

function index_socr(model::Model, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	ibmdp = index_bmdp(model, env, quat = quat)
	ibfri = index_bfri(model, env, quat = quat)
	isocr = [ibmdp; ibfri]
	return isocr
end

function index_socri(model::Model, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
	isocri = Vector{Vector{Int}}()
	return isocri
end

function index_socri(model::Model, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	isocri = Vector{Vector{Int}}()
	return isocri
end

# struct DynamicsStructure <: IndicesStructure 
# 	ny::Int 
# 	dynz::Vector{Int} 
# 	dynΔ::Vector{Int} 
# 	dyn::Vector{Int} 
# 	rst::Vector{Int} 
# 	alt::Vector{Int} 
# end

function IndicesOptimization(model::Model, env::Environment)
	# Dimensions
	nz = num_var(model, env, quat = false)
	nΔ = num_var(model, env, quat = true)
	ny = num_bilinear(model, env)

	# Variables
	dynz = index_q2(model, env, quat = false)
	dynΔ = index_q2(model, env, quat = true)
	ortz = index_ort(model, env, quat = false)
	ortΔ = index_ort(model, env, quat = true)
	socz = index_soc(model, env, quat = false)
	socΔ = index_soc(model, env, quat = true)

	# Residual
	# dyn = [dyn]
	# rst = [s1  - ..., ≡ ialt
	#        η1  - ...,
	#        s2  - ...,]
	# bil = [γ1 .* s1 .- κ;
	#        b1 .* η1 .- κ;
	#        ψ1 .* s2 .- κ]

	equr = index_equr(model, env, quat = true)
	ortr = index_ortr(model, env, quat = true)
	socr = index_socr(model, env, quat = true)
	socri = index_socri(model, env, quat = true)

	dyn = index_dyn(model, env, quat = true)
	imp = index_imp(model, env, quat = true)
	mdp = index_mdp(model, env, quat = true)
	fri = index_fri(model, env, quat = true)
	bimp = index_bimp(model, env, quat = true)
	bmdp = index_bmdp(model, env, quat = true)
	bfri = index_bfri(model, env, quat = true)

	rst = [imp; mdp; fri]
	bil = [bimp; bmdp; bfri]
<<<<<<< HEAD
	alt = imp

	s = OptimizationIndices(
		nz, nΔ, ny, nquat,
		dynz, dynΔ, ortz, ortΔ, socz, socΔ,
		equr, ortr, socr, dyn, rst, bil, alt)
	return s
end
=======
	
	IndicesOptimization(
		nz, nΔ,
		ortz, ortΔ, socz, socΔ,
		equr, ortr, socr, socri, bil)
end
>>>>>>> be9e65f... simulator from RoboDojo
