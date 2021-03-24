"""
	ImplicitTraj{T}
This structure holds the trajectory of evaluations and Jacobians of the implicit dynamics.
These evaluations and Jacobians are computed using a linearizedimation computed around `lin`.
"""
struct ImplicitTraj{T}
	lin::Vector{LinearizedStep{T}}
	d::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}} # dynamics violation
	dq2::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	dγ1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	db1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	δq0::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # q0 solution gradient length=H
	δq1::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # q1 solution gradient length=H
	δu1::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # u1 solution gradient length=H
	ip::Vector{InteriorPoint{T}}
end

function ImplicitTraj(ref_traj::ContactTraj, model::ContactDynamicsModel;
	κ = ref_traj.κ[1],
	opts = InteriorPointOptions(
			κ_init = κ[1],
			κ_tol = 2.0 * κ[1],
			r_tol = 1.0e-8,
			diff_sol = true))

	H = ref_traj.H

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq + nc + nb
	nz = num_var(model)
	nθ = num_data(model)

	lin = [LinearizedStep(model, ref_traj.z[t], ref_traj.θ[t], κ) for t = 1:H]

	ip =  [interior_point(zeros(num_var(model)), zeros(num_data(model)),
			 idx_ineq = inequality_indices(model),
			 r! = lin[t].methods.r!,
			 rz! = lin[t].methods.rz!,
			 rθ! = lin[t].methods.rθ!,
			 rz = model.spa.rz_sp,
			 rθ = model.spa.rθ_sp,
			 opts = opts) for t = 1:H]

	# views
	d = [view(ip[t].z, 1:nd) for t = 1:H]
	dq2 = [view(ip[t].z, 1:nq) for t = 1:H]
	dγ1 = [view(ip[t].z, nq .+ (1:nc)) for t = 1:H]
	db1 = [view(ip[t].z, nq + nc .+ (1:nb)) for t = 1:H]

	off = 0
	δq0 = [view(ip[t].δz, 1:nd, off .+ (1:nq)) for t = 1:H]; off += nq
	δq1 = [view(ip[t].δz, 1:nd, off .+ (1:nq)) for t = 1:H]; off += nq
	δu1 = [view(ip[t].δz, 1:nd, off .+ (1:nu)) for t = 1:H]; off += nu

	return ImplicitTraj(lin, d, dq2, dγ1, db1, δq0, δq1, δu1, ip)
end

"""
	implicit_dynamics!(im_traj::ImplicitTraj, model::ContactDynamicsModel,
		traj::ContactTraj; κ = traj.κ)
Compute the evaluations and Jacobians of the implicit dynamics on the trajectory 'traj'. The computation is
linearizedimated since it relies on a linearization about a reference trajectory.
"""
function implicit_dynamics!(im_traj::ImplicitTraj, model::ContactDynamicsModel,
	traj::ContactTraj; κ = traj.κ) where {T, D}

	for t = 1:traj.H
		#@assert abs.(im_traj.lin[t].κ - κ[1]) / κ[1] < 1.0e-5 # check that the κ are consistent between the optimized trajectory (traj)
		# and the linearization (impl.lin).

		# initialized solver
		z_initialize!(im_traj.ip[t].z, model, traj.q[t+2]) #TODO: try alt. schemes
		im_traj.ip[t].θ .= traj.θ[t]

		# solve
		status = interior_point!(im_traj.ip[t])
		!status && (@warn "implicit dynamics failure (t = $t)")

		# compute dynamics violation
		im_traj.dq2[t] .-= traj.q[t+2]
		im_traj.dγ1[t] .-= traj.γ[t]
		im_traj.db1[t] .-= traj.b[t]
	end

	return nothing
end
