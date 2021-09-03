# TODO
# we need to have nquat directly inside the model


"""
	Returns the indices of q2 in z and Δq2 in Δz. Important to note that for the indexing
	in z, we want to set nquat = 0. For indexing in Δz, we set nquat to the
	number of quaternions.
"""
function index_q2(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	iq2 = Vector(1:nq - nquat)
	return iq2
end

"""
	Returns the indices of γ1 in z and Δγ1 in Δz. Important to note that for the indexing
	in z, we want to set nquat = 0. For indexing in Δz, we set nquat to the
	number of quaternions.
"""
function index_γ1(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	off = (nq - nquat)
	iγ1 = Vector(off .+ (1:nc))
	return iγ1
end

"""
	Returns the indices of b1 in z and Δb1 in Δz. Important to note that for the indexing
	in z, we want to set nquat = 0. For indexing in Δz, we set nquat to the
	number of quaternions.
"""
function index_b1(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc
	ib1 = Vector(off .+ (1:nb))
	return ib1
end

"""
	Returns the indices of ψ1 in z and Δψ1 in Δz. Important to note that for the indexing
	in z, we want to set nquat = 0. For indexing in Δz, we set nquat to the
	number of quaternions.
"""
function index_ψ1(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb
	iψ1 = Vector(off .+ (1:nc))
	return iψ1
end

"""
	Returns the indices of s1 in z and Δs1 in Δz. Important to note that for the indexing
	in z, we want to set nquat = 0. For indexing in Δz, we set nquat to the
	number of quaternions.
"""
function index_s1(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc
	is1 = Vector(off .+ (1:nc))
	return is1
end

"""
	Returns the indices of η1 in z and Δη1 in Δz. Important to note that for the indexing
	in z, we want to set nquat = 0. For indexing in Δz, we set nquat to the
	number of quaternions.
"""
function index_η1(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc + nc
	iη1 = Vector(off .+ (1:nb))
	return iη1
end

"""
	Returns the indices of s2 in z and Δs2 in Δz. Important to note that for the indexing
	in z, we want to set nquat = 0. For indexing in Δz, we set nquat to the
	number of quaternions.
"""
function index_s2(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)
	off = (nq - nquat) + nc + nb + nc + nc + nb
	is2 = Vector(off .+ (1:nc))
	return is2
end


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
function linearization_var_index(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)

	ny = nc + nb + nc

	off = 0
	iw1 = off .+ Vector(1:nq - nquat); off += nq - nquat
	iw2 = off .+ Vector(1:ny); off += ny
	iw3 = off .+ Vector(1:ny); off += ny
	return iw1, iw2, iw3
end

"""
	Returns the indices of the residual r, where 3 groups are formed.
	r = (rdyn, rrst, rbil)
"""
function linearization_term_index(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)

	ny = 2nc + nb
	# dyn = [dyn]
	# rst = [s1  - ..., ≡ ialt
	#        η1  - ...,
	#        s2  - ...,]
	# bil = [γ1 .* s1 .- κ;
	#        b1 .* η1 .- κ;
	#        ψ1 .* s2 .- κ]
	idyn  = Vector(1:nq - nquat)
	irst1 = Vector(nq - nquat .+ (1:nc))
	irst2 = Vector(nq - nquat + nc .+ (1:nb))
	irst3 = Vector(nq - nquat + nc + nb .+ (1:nc))
	irst  = [irst1; irst2; irst3]
	ialt  = irst1
	ibil1 = Vector(nq - nquat + nb + 2nc .+ (1:nc))
	ibil2 = Vector(nq - nquat + nb + 3nc .+ (1:nb))
	ibil3 = Vector(nq - nquat + nb + 3nc + nb .+ (1:nc))
	ibil  = [ibil1; ibil2; ibil3]
	return idyn, irst, ibil, ialt
end

"""
	Returns the 3 residual bilinear terms and the 3 couples of variables
	associated with them. For each bilinear residual terms we associate 2
	variables.
"""
function get_bilinear_indices(model::ContactModel, env::Environment; nquat::Int = 0)
	nq = model.dim.q
	nc = model.dim.c
	nb = nc * friction_dim(env)

	terms = [SVector{nc,Int}(nq - nquat + nb + 2nc .+ (1:nc)),      # γ1, s1
			 SVector{nb,Int}(nq - nquat + nb + 3nc .+ (1:nb)),      # b1, η1
			 SVector{nc,Int}(nq - nquat + nb + 3nc + nb .+ (1:nc))] # ψ1, s2

	vars = [[SVector{nc,Int}(nq - nquat .+ (1:nc)),
			 SVector{nc}(nq - nquat + 2nc +  nb .+ (1:nc))], # γ1, s1
			[SVector{nb,Int}(nq - nquat + nc .+ (1:nb)),
			 SVector{nb}(nq - nquat + 3nc +  nb .+ (1:nb))], # b1, η1
			[SVector{nc,Int}(nq - nquat + nc + nb .+ (1:nc)),
			 SVector{nc}(nq - nquat + 3nc + 2nb .+ (1:nc))], # ψ1, s2
			]
	return terms, vars
end

"""
	Returns the positive orthant indices in z or Δz.
"""
function inequality_indices(model::ContactModel, env::Environment{<:World,LinearizedCone}; nquat::Int = 0)
	iγ1 = index_γ1(model, env, nquat = nquat)
	ib1 = index_b1(model, env, nquat = nquat)
	iψ1 = index_ψ1(model, env, nquat = nquat)
	is1 = index_s1(model, env, nquat = nquat)
	iη1 = index_η1(model, env, nquat = nquat)
	is2 = index_s2(model, env, nquat = nquat)
	return [iγ1; ib1; iψ1; is1; iη1; is2]
end

"""
	Returns the positive orthant indices in z or Δz.
"""
function inequality_indices(model::ContactModel, env::Environment{<:World,NonlinearCone}; nquat::Int = 0)
	# nb = model.dim.c * friction_dim(env)
	# collect([(model.dim.q .+ (1:model.dim.c))...,
	#          (model.dim.q + model.dim.c + nb + nb + model.dim.c .+ (1:model.dim.c))...])
	iγ1 = index_γ1(model, env, nquat = nquat)
	is1 = index_s1(model, env, nquat = nquat)
	return [iγ1; is1]
end

"""
	Returns the second order cone indices in z or Δz.
"""
soc_indices(model::ContactModel, env::Environment{<:World,LinearizedCone}; nquat::Int = 0) = Vector{Int}[]

"""
	Returns the second order cone indices in z or Δz.
"""
function soc_indices(model::ContactModel, env::Environment{<:World,NonlinearCone}; nquat::Int = 0)
	nc = model.dim.c
	nf = friction_dim(env)

	ib1 = index_b1(model, env, nquat = nquat) # primal cones: vector part
	iψ1 = index_ψ1(model, env, nquat = nquat) # primal cones: scalar part
	iη1 = index_η1(model, env, nquat = nquat) # dual cones: vector part
	is2 = index_s2(model, env, nquat = nquat) # dual cones: scalar part

	# b_idx = nq + nc .+ (1:nb)
	# η_idx = nq + nc + nb .+ (1:(nb + nc))
	# s2_idx = nq + nc + nb + nb + nc + nc .+ (1:nc)
	#
	# pr_idx = [[s2_idx[i]; b_idx[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
	# du_idx = [[η_idx[(i - 1) * ne .+ (1:ne)]...] for i = 1:nc]
	#
	# [pr_idx..., du_idx...]
	pr_idx = [[iψ1[i]; ib1[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
	du_idx = [[is2[i]; iη1[(i - 1) * nf .+ (1:nf)]] for i = 1:nc]
	[pr_idx..., du_idx...]
end
