"""
	ImplicitTraj{T,nd}
This structure holds the trajectory of evaluations and Jacobians of the implicit dynamics.
These evaluations and Jacobians are computed using a linear approximation computed around `lin`.

"""
mutable struct ImplicitTraj{T,nd}
	H::Int                               # horizon length
	lin::Vector{LinStep{T}}              # linearization point  length=H
	d::Vector{SizedVector{nd,T}}         # dynamics violation   length=H
	δz::Vector{SparseMatrixCSC{T,Int}}   # solution gradient    length=H
	δq0::Vector{SparseMatrixCSC{T,Int}}  # q0 solution gradient length=H
	δq1::Vector{SparseMatrixCSC{T,Int}}  # q0 solution gradient length=H
	δu1::Vector{SparseMatrixCSC{T,Int}}  # q0 solution gradient length=H
end

"""
	Create dummy ImplicitTraj.
"""
function ImplicitTraj(H::Int, model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq+nc+nb
	nz = num_var(model)
	nθ = num_data(model)
	lin = [LinStep(model) for k = 1:H]
	d = [zeros(SizedVector{nd}) for k=1:H]
	δz  = [spzeros(nz,nθ) for k=1:H]
	δq0 = [spzeros(nd,nq) for k=1:H]
	δq1 = [spzeros(nd,nq) for k=1:H]
	δu1 = [spzeros(nd,nu) for k=1:H]
	return ImplicitTraj{eltype(d[1]),nd}(H,lin,d,δz,δq0,δq1,δu1)
end

"""
	linearization!(model::ContactDynamicsModel, ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ},
		impl::ImplicitTraj{T,nd}, κ::T=ref_traj.κ) where {T,nq,nu,nw,nc,nb,nz,nθ,nd}
Linearization of the model around the reference trajectory `ref_traj`, the resulting linearization is stored
in `impl.lin`.
"""
function linearization!(model::ContactDynamicsModel, ref_traj,
	impl::ImplicitTraj{T, nd}, κ=ref_traj.κ) where {T,nd}
	@assert impl.H == ref_traj.H
	H = ref_traj.H
	for k = 1:H
		# Compute linearization
		# we need to figure out, the correct way to handle ref_traj that are optimized at κ=1e-10,
		# and linearization at κ=1e-3
		impl.lin[k] = LinStep(model, ref_traj.z[k], ref_traj.θ[k], κ[1])
	end
	return nothing
end

"""
		implicit_dynamics!(model::ContactDynamicsModel, traj::ContactTraj{T,nq,nu,nw,nc,nb},
	impl::ImplicitTraj{T,nd}; κ::T=traj.κ) where {T,nq,nu,nw,nc,nb,nd}
Compute the evaluations and Jacobians of the implicit dynamics on the trajectory 'traj'. The computation is
approximated since it relies on a linearization about a reference trajectory.
"""
function implicit_dynamics!(model::ContactDynamicsModel, traj::ContactTraj,
	impl::ImplicitTraj{T, nd}; κ=traj.κ) where {T,nd}
	@assert impl.H >= traj.H
	H = traj.H
	# Compute the implicit dynamics
	# constraint violation (solve linearized contact problem)
	# constraint gradient (implicit function theorem)

	# INITIALIZATION
	ip = interior_point(zeros(num_var(model)), zeros(num_data(model)),
		idx_ineq = inequality_indices(model),
		r! = model.res.r,
		rz! = model.res.rz,
		rθ! = model.res.rθ,
		rz = model.spa.rz_sp,
		rθ = model.spa.rθ_sp,
		opts = InteriorPointOptions(
			κ_init=κ[1],
			κ_tol=2κ[1],
			r_tol=1e-8,
			diff_sol=true)
		)

	for k = 1:H
		@assert abs.(impl.lin[k].κ0 - κ[1])/κ[1] < 1e-5 # check that the κ are consistent between the optimized trajectory (traj)
		# and the linearization (impl.lin).

		# Define local residual functions
		# residual
		r!(r, z, θ, κ) = model.approx.r(r, z, θ, κ, impl.lin[k].z0, impl.lin[k].θ0, impl.lin[k].r0, impl.lin[k].rz0, impl.lin[k].rθ0)
		# residual Jacobian wrt z
		rz!(rz, z, θ) = model.approx.rz(rz, z, impl.lin[k].rz0)
		# residual Jacobian wrt θ
		rθ!(rθ, z, θ) = model.approx.rθ(rθ, impl.lin[k].rθ0)
		# Set the residual functions
		ip.methods.r! = r!
		ip.methods.rz! = rz!
		ip.methods.rθ! = rθ!

		z_initialize!(ip.z, model, traj.q[k+2]) # initialize with our best guess.
		#################################
		# ip.z .= copy(traj.z[k]) # Maybe better needs testing
		#################################
		ip.θ .= traj.θ[k]
		status = interior_point!(ip)

		q2 = traj.q[k+2]
		γ1 = traj.γ[k]
		b1 = traj.b[k]

		nq = model.dim.q
		nu = model.dim.u

		impl.d[k] = ip.z[1:nd] - [q2; γ1; b1]
		impl.δz[k]  = copy(ip.δz)
		off = 0
		# TODO we need to have a struct that stores the indices of q2, γ1, b1 in z
		# TODO we need to have a struct that stores the indices of q0, q1, u1 in θ
		impl.δq0[k] = copy(ip.δz[1:nd, off .+ (1:nq)]); off += nq # ndxnq
		impl.δq1[k] = copy(ip.δz[1:nd, off .+ (1:nq)]); off += nq # ndxnq
		impl.δu1[k] = copy(ip.δz[1:nd, off .+ (1:nu)]); off += nu # ndxnu
	end
	return nothing
end
