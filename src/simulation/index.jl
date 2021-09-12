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
function index_q2(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nquat = quat ? model.dim.quat : 0
	iq2 = Vector(1:nq - nquat)
	return iq2
end

"""
	Returns the indices of γ1 in z and Δγ1 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_γ1(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
	off = (nq - nquat)
	iγ1 = Vector(off .+ (1:nc))
	return iγ1
end

"""
	Returns the indices of b1 in z and Δb1 in Δz. Important to note that for the indexing
	in z, we want to set quat = false. For indexing in Δz, we set quat = true, so that nquat =
	number of quaternions.
"""
function index_b1(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
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
function index_ψ1(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
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
function index_s1(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
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
function index_η1(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
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
function index_s2(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
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
function index_q0(model::ContactModel)
	nq = model.dim.q
	iq0 = Vector(1:nq)
	return iq0
end

"""
	Returns the indices of q1 in θ.
"""
function index_q1(model::ContactModel)
	nq = model.dim.q
	off = nq
	iq1 = Vector(off .+ (1:nq))
	return iq1
end

"""
	Returns the indices of u1 in θ.
"""
function index_u1(model::ContactModel)
	nq = model.dim.q
	nu = model.dim.u
	off = nq + nq
	iu1 = Vector(off .+ (1:nu))
	return iu1
end

"""
	Returns the indices of w1 in θ.
"""
function index_w1(model::ContactModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	off = nq + nq + nu
	iw1 = Vector(off .+ (1:nw))
	return iw1
end

"""
	Returns the indices of μ in θ.
"""
function index_μ(model::ContactModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	off = nq + nq + nu + nw
	iμ = Vector(off .+ (1:1))
	return iμ
end

"""
	Returns the indices of h in θ.
"""
function index_h(model::ContactModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
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
function index_dyn(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nquat = quat ? model.dim.quat : 0
	idyn = Vector(1:nq - nquat)
	return idyn
end

"""
	Returns the indices of the impact equation in the residual r.
"""
function index_imp(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
	off = (nq - nquat)
	iimp = Vector(off .+ (1:nc))
	return iimp
end

"""
	Returns the indices of the maximum dissipation principle equation in the residual r.
"""
function index_mdp(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc
	imdp = Vector(off .+ (1:nb))
	return imdp
end

"""
	Returns the indices of the friction cone equation in the residual r.
"""
function index_fri(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb
	ifri = Vector(off .+ (1:nc))
	return ifri
end

"""
	Returns the indices of the bilinear impact equation in the residual r.
"""
function index_bimp(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc
	ibimp = Vector(off .+ (1:nc))
	return ibimp
end

"""
	Returns the indices of the bilinear maximum dissipation principle equation in the residual r.
"""
function index_bmdp(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc + nc
	ibmdp = Vector(off .+ (1:nb))
	return ibmdp
end

"""
	Returns the indices of the bilinear friction equation in the residual r.
"""
function index_bfri(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
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
function linearization_var_index(model::ContactModel, env::Environment; quat::Bool = false)
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
function linearization_term_index(model::ContactModel, env::Environment; quat::Bool = false)
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
	Returns the 3 residual bilinear terms and the 3 couples of variables
	associated with them. For each bilinear residual terms we associate 2
	variables.
"""
function get_bilinear_indices(model::ContactModel, env::Environment; quat::Bool = false)
	idyn = index_dyn(model, env, quat = quat)
	iimp = index_imp(model, env, quat = quat)
	imdp = index_mdp(model, env, quat = quat)
	ifri = index_fri(model, env, quat = quat)
	ibimp = index_bimp(model, env, quat = quat)
	ibmdp = index_bmdp(model, env, quat = quat)
	ibfri = index_bfri(model, env, quat = quat)

	iq2 = index_q2(model, env, quat = quat)
	iγ1 = index_γ1(model, env, quat = quat)
	ib1 = index_b1(model, env, quat = quat)
	iψ1 = index_ψ1(model, env, quat = quat)
	is1 = index_s1(model, env, quat = quat)
	iη1 = index_η1(model, env, quat = quat)
	is2 = index_s2(model, env, quat = quat)

	terms = [SVector{length(ibimp),Int}(ibimp), # γ1, s1
			 SVector{length(ibmdp),Int}(ibmdp), # b1, η1
			 SVector{length(ibfri),Int}(ibfri)] # ψ1, s2

	vars = [[SVector{length(iγ1),Int}(iγ1),
			 SVector{length(is1),Int}(is1)], # γ1, s1
			[SVector{length(ib1),Int}(ib1),
			 SVector{length(iη1),Int}(iη1)], # b1, η1
			[SVector{length(iψ1),Int}(iψ1),
			 SVector{length(is2),Int}(is2)], # ψ1, s2
			]
	return terms, vars
end

"""
	Returns the positive orthant indices in z or Δz.
"""
function inequality_indices(model::ContactModel, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
	iγ1 = index_γ1(model, env, quat = quat)
	ib1 = index_b1(model, env, quat = quat)
	iψ1 = index_ψ1(model, env, quat = quat)
	is1 = index_s1(model, env, quat = quat)
	iη1 = index_η1(model, env, quat = quat)
	is2 = index_s2(model, env, quat = quat)
	return [iγ1; ib1; iψ1; is1; iη1; is2]
end

"""
	Returns the positive orthant indices in z or Δz.
"""
function inequality_indices(model::ContactModel, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	iγ1 = index_γ1(model, env, quat = quat)
	is1 = index_s1(model, env, quat = quat)
	return [iγ1; is1]
end

"""
	Returns the positive orthant indices in z or Δz.
"""
function index_ort(model::ContactModel, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
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
function index_ort(model::ContactModel, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	iγ1 = index_γ1(model, env, quat = quat)
	is1 = index_s1(model, env, quat = quat)
	return [iγ1, is1]
end

"""
	Returns the second order cone indices in z or Δz.
"""
index_soc(model::ContactModel, env::Environment{<:World,LinearizedCone}; quat::Bool = false) = [Vector{Int}[], Vector{Int}[]]

"""
	Returns the second order cone indices in z or Δz.
"""
function index_soc(model::ContactModel, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	nc = model.dim.c
	nf = friction_dim(env)

	ib1 = index_b1(model, env, quat = quat) # primal cones: vector part
	iψ1 = index_ψ1(model, env, quat = quat) # primal cones: scalar part
	iη1 = index_η1(model, env, quat = quat) # dual cones: vector part
	is2 = index_s2(model, env, quat = quat) # dual cones: scalar part

	pr_idx = [[iψ1[i]; ib1[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
	du_idx = [[is2[i]; iη1[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
	[pr_idx, du_idx]
end

function num_var(model::ContactModel, env::Environment; quat::Bool = false)
	nq = model.dim.q
	nc = model.dim.c
	nquat = quat ? model.dim.quat : 0
	nb = nc * friction_dim(env)
	(nq - nquat) + nc + nb + nc + nc + nb + nc
end

function num_data(model::ContactModel)
	nq = model.dim.q
	nu = model.dim.u
	nw = model.dim.w
	nq + nq + nu + nw + 1 + 1
end

function num_bilinear(model::ContactModel, env::Environment)
	nc = model.dim.c
	nb = nc * friction_dim(env)
	nc + nb + nc
end

################################################################################
# Packing and Unpacking Variables
################################################################################

function unpack_θ(model::ContactModel, θ)
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

function pack_θ(model::ContactModel, q0, q1, u1, w1, μ, h)
	return [q0; q1; u1; w1; μ; h]
end

# function pack_θ(model::ContactModel, q0, q1, u1, w1, μ, h)
# 	nθ = num_data(model)
#
# 	iq0 = index_q0(model)
# 	iq1 = index_q1(model)
# 	iu1 = index_u1(model)
# 	iw1 = index_w1(model)
# 	iμ  = index_μ(model)
# 	ih  = index_h(model)
#
# 	θ = zeros(nθ)
# 	θ[iq0] = q0
# 	θ[iq1] = q1
# 	θ[iu1] = u1
# 	θ[iw1] = w1
# 	θ[iμ] = μ
# 	θ[ih] = h
# 	return θ
# end

function unpack_z(model::ContactModel, env::Environment, z)
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

function pack_z(model::ContactModel, env::Environment{<:World,LinearizedCone}, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_func(model, env, q2)
	s2 = model.μ_world * γ1 .- E_func(model, env) * b1
	return pack_z(model, env, q2, γ1, b1, ψ1, s1, η1, s2)
end

function pack_z(model::ContactModel, env::Environment{<:World,NonlinearCone}, q2, γ1, b1, ψ1, η1)
	s1 = ϕ_func(model, q2)
	s2 = model.μ_world .* γ1
	return pack_z(model, env, q2, γ1, b1, ψ1, s1, η1, s2)
end

function pack_z(model::ContactModel, env::Environment, q2, γ1, b1, ψ1, s1, η1, s2)
	return [q2; γ1; b1; ψ1; s1; η1; s2]
end


# function pack_z(model::ContactModel, env::Environment, q2, γ1, b1, ψ1, s1, η1, s2)
# 	nz = num_var(model, env)
# 	iq2 = index_q2(model, env, quat = false)
# 	iγ1 = index_γ1(model, env, quat = false)
# 	ib1 = index_b1(model, env, quat = false)
# 	iψ1 = index_ψ1(model, env, quat = false)
# 	is1 = index_s1(model, env, quat = false)
# 	iη1 = index_η1(model, env, quat = false)
# 	is2 = index_s2(model, env, quat = false)
#
# 	z = zeros(nz)
#
# 	z[iq2] = q2
# 	z[iγ1] = γ1
# 	z[ib1] = b1
# 	z[iψ1] = ψ1
# 	z[is1] = s1
# 	z[iη1] = η1
# 	z[is2] = s2
# 	return z
# end


################################################################################
# Aggregated indices
################################################################################

function index_eq(model::ContactModel, env::Environment; quat::Bool = false)
	iq2 = index_q2(model, env, quat = quat)
	iγ1 = index_γ1(model, env, quat = quat)
	ib1 = index_b1(model, env, quat = quat)
	iψ1 = index_ψ1(model, env, quat = quat)

	iequ = [iq2; iγ1; ib1; iψ1]
	iort = index_ort(model, env, quat = quat)
	isoc = index_soc(model, env, quat = quat)
	return iequ, iort, isoc
end

function index_variable(model::ContactModel, env::Environment; quat::Bool = false)
	iq2 = index_q2(model, env, quat = quat)
	iγ1 = index_γ1(model, env, quat = quat)
	ib1 = index_b1(model, env, quat = quat)
	iψ1 = index_ψ1(model, env, quat = quat)

	iequ = [iq2; iγ1; ib1; iψ1]
	iort = index_ort(model, env, quat = quat)
	isoc = index_soc(model, env, quat = quat)
	return iequ, iort, isoc
end

function index_residual(model::ContactModel, env::Environment; quat::Bool = false)
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
	ibil_ort = index_ortr(model, env, quat = quat)
	ibil_soc = index_socr(model, env, quat = quat)
	return idyn, irst, ibil, ialt, ibil_ort, ibil_soc
end

function index_equr(model::ContactModel, env::Environment; quat::Bool = false)
	idyn = index_dyn(model, env, quat = quat)
	iimp = index_imp(model, env, quat = quat)
	imdp = index_mdp(model, env, quat = quat)
	ifri = index_fri(model, env, quat = quat)
	iequ = [idyn; iimp; imdp; ifri]
	return iequ
end

function index_ortr(model::ContactModel, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
	ibimp = index_bimp(model, env, quat = quat)
	ibmdp = index_bmdp(model, env, quat = quat)
	ibfri = index_bfri(model, env, quat = quat)
	iortr = [ibimp; ibmdp; ibfri]
	return iortr
end

function index_ortr(model::ContactModel, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	ibimp = index_bimp(model, env, quat = quat)
	iortr = ibimp
	return iortr
end

function index_socr(model::ContactModel, env::Environment{<:World,LinearizedCone}; quat::Bool = false)
	isocr = Vector{Int64}()
	return isocr
end

function index_socr(model::ContactModel, env::Environment{<:World,NonlinearCone}; quat::Bool = false)
	ibmdp = index_bmdp(model, env, quat = quat)
	ibfri = index_bfri(model, env, quat = quat)
	isocr = [ibmdp; ibfri]
	return isocr
end

################################################################################
# Optimization Space
################################################################################

abstract type OptimizationSpace
end

mutable struct OptimizationSpace12 <: OptimizationSpace
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
	socr::Vector{Vector{Int}} # indices of the residual associated with the Second Order Cone constraints in r
	dyn::Vector{Int} # indices of the residual associated with the dynamics constraints in r
	rst::Vector{Int} # indices of the residual associated with the remaining constraints in r
	bil::Vector{Int} # indices of the residual associated with the bilinear constraints in r
	alt::Vector{Int} # indices of the residual associated with the altitude constraints in r
end

function OptimizationSpace12(model::ContactModel, env::Environment)
	# Dimensions
	nz = num_var(model, env, quat = false)
	nΔ = num_var(model, env, quat = true)
	ny = num_bilinear(model, env)
	nquat = model.dim.quat

	# Variables
	dynz = index_q2(model, env, quat = false)
	dynΔ = index_q2(model, env, quat = true)
	ortz = index_ort(model, env, quat = true)
	ortΔ = index_ort(model, env, quat = false)
	socz = index_soc(model, env, quat = true)
	socΔ = index_soc(model, env, quat = false)

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

	dyn = index_dyn(model, env, quat = true)
	imp = index_imp(model, env, quat = true)
	mdp = index_mdp(model, env, quat = true)
	fri = index_fri(model, env, quat = true)
	bimp = index_bimp(model, env, quat = true)
	bmdp = index_bmdp(model, env, quat = true)
	bfri = index_bfri(model, env, quat = true)

	rst = [imp; mdp; fri]
	bil = [bimp; bmdp; bfri]
	alt = imp

	s = OptimizationSpace12(
		nz, nΔ, ny, nquat,
		dynz, dynΔ, ortz, ortΔ, socz, socΔ,
		equr, ortr, socr, dyn, rst, bil, alt)
	return s
end

function OptimizationSpace12()
	v1 = Vector{Int}()
	v2 = Vector{Vector{Int}}()
	v3 = Vector{Vector{Vector{Int}}}()

	s = OptimizationSpace12(
		0, 0, 0, 0,
		v1, v1, v2, v2, v3, v3,
		v1, v1, v2, v1, v1, v1, v1)
	return s
end
