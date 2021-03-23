"""
	ImplicitTraj{T,D}
This structure holds the trajectory of evaluations and Jacobians of the implicit dynamics.
These evaluations and Jacobians are computed using a linearizedimation computed around `lin`.

"""
struct ImplicitTraj{T,D}
	lin::Vector{LinearizedStep{T}}       # linearization point  length=H
	d::Vector{SizedArray{Tuple{D},T,1,1,Array{T,1}}} # dynamics violation   length=H
	δz::Vector{SparseMatrixCSC{T,Int}}   # solution gradient    length=H
	δq0::Vector{SparseMatrixCSC{T,Int}}  # q0 solution gradient length=H
	δq1::Vector{SparseMatrixCSC{T,Int}}  # q0 solution gradient length=H
	δu1::Vector{SparseMatrixCSC{T,Int}}  # q0 solution gradient length=H
end

function implicit_trajectory(H::Int, model::ContactDynamicsModel)
	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq + nc + nb
	nz = num_var(model)
	nθ = num_data(model)

	lin = [LinearizedStep(model) for k = 1:H]
	d = [zeros(SizedVector{nd}) for k = 1:H]
	δz  = [spzeros(nz, nθ) for k = 1:H]
	δq0 = [spzeros(nd, nq) for k = 1:H]
	δq1 = [spzeros(nd, nq) for k = 1:H]
	δu1 = [spzeros(nd, nu) for k = 1:H]

	return ImplicitTraj(lin, d, δz, δq0, δq1, δu1)
end

"""
	linearization!(im_traj::ImplicitTraj, model::ContactDynamicsModel,
		ref_traj::ContactTraj; κ = ref_traj.κ)
Linearization of the model around the reference trajectory `ref_traj`, the resulting linearization is stored
in `im_traj.lin`.
"""
function linearization!(im_traj::ImplicitTraj, model::ContactDynamicsModel,
	ref_traj::ContactTraj; κ = ref_traj.κ)

	H = ref_traj.H

	for t = 1:H
		# Compute linearization
		im_traj.lin[t] = LinearizedStep(model, ref_traj.z[t], ref_traj.θ[t], κ[1])
	end

	return nothing
end

"""
	implicit_dynamics!(im_traj::ImplicitTraj, model::ContactDynamicsModel,
		traj::ContactTraj; κ = traj.κ)
Compute the evaluations and Jacobians of the implicit dynamics on the trajectory 'traj'. The computation is
linearizedimated since it relies on a linearization about a reference trajectory.
"""
function implicit_dynamics!(im_traj::ImplicitTraj{T, D}, model::ContactDynamicsModel,
	traj::ContactTraj; κ = traj.κ) where {T, D}

	H = traj.H
	# Compute the implicit dynamics

	# INITIALIZATION
	ip = interior_point(zeros(num_var(model)), zeros(num_data(model)),
		idx_ineq = inequality_indices(model),
		r! = model.res.r!,
		rz! = model.res.rz!,
		rθ! = model.res.rθ!,
		rz = model.spa.rz_sp,
		rθ = model.spa.rθ_sp,
		opts = InteriorPointOptions(
			κ_init = κ[1],
			κ_tol = 2κ[1],
			r_tol = 1e-8,
			diff_sol = true))

	for t = 1:H
		@assert abs.(im_traj.lin[t].κ - κ[1]) / κ[1] < 1e-5 # check that the κ are consistent between the optimized trajectory (traj)
		# and the linearization (impl.lin).

		# Set the residual functions
		ip.methods.r! = im_traj.lin[t].methods.r!
		ip.methods.rz! = im_traj.lin[t].methods.rz!
		ip.methods.rθ! = im_traj.lin[t].methods.rθ!

		z_initialize!(ip.z, model, traj.q[t+2]) # initialize with our best guess.
		#################################
		# ip.z .= copy(traj.z[k]) # Maybe better needs testing
		#################################
		ip.θ .= traj.θ[t]

		status = interior_point!(ip)

		q2 = traj.q[t+2]
		γ1 = traj.γ[t]
		b1 = traj.b[t]

		nq = model.dim.q
		nu = model.dim.u

		im_traj.d[t] = ip.z[1:D] - [q2; γ1; b1] # TODO fix
		im_traj.δz[t]  = copy(ip.δz)
		off = 0
		# TODO we need to have a struct that stores the indices of q2, γ1, b1 in z
		# TODO we need to have a struct that stores the indices of q0, q1, u1 in θ
		im_traj.δq0[t] = copy(ip.δz[1:D, off .+ (1:nq)]); off += nq # ndxnq
		im_traj.δq1[t] = copy(ip.δz[1:D, off .+ (1:nq)]); off += nq # ndxnq
		im_traj.δu1[t] = copy(ip.δz[1:D, off .+ (1:nu)]); off += nu # ndxnu
	end

	return nothing
end
