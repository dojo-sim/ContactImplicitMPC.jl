"""
	ImplicitTraj{T}
This structure holds the trajectory of evaluations and Jacobians of the implicit dynamics.
These evaluations and Jacobians are computed using a linearizedimation computed around `lin`.
"""
mutable struct ImplicitTraj{T}
	lin::Vector{LinearizedStep{T}}
	d::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}} # dynamics violation
	dq2::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	dγ1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	db1::Vector{SubArray{T,1,Array{T,1},Tuple{UnitRange{Int}},true}}
	# δq0::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # q0 solution gradient length=H
	# δq1::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # q1 solution gradient length=H
	# δu1::Vector{SubArray{T,2,SparseArrays.SparseMatrixCSC{T,Int},Tuple{UnitRange{Int},UnitRange{Int}},false}}  # u1 solution gradient length=H
	δq0::Vector{SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}} # q0 solution gradient length=H
	δq1::Vector{SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}}  # q1 solution gradient length=H
	δu1::Vector{SubArray{Float64,2,Array{Float64,2},Tuple{UnitRange{Int64},UnitRange{Int64}},false}}  # u1 solution gradient length=H
	ip::Vector{InteriorPoint{T}}
end


function ImplicitTraj(ref_traj::ContactTraj, model::ContactDynamicsModel;
	κ = ref_traj.κ[1],
	max_time = 60.0,
	opts = InteriorPointOptions(
			κ_init = κ[1],
			κ_tol = 2.0 * κ[1],
			r_tol = 1.0e-8,
			diff_sol = true,
			solver=:empty_solver,
			max_time=max_time))

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
			 r! = r!,
			 rz! = rz!,
			 rθ! = rθ!,
			 r  = RLin(model, lin[t].z, lin[t].θ, lin[t].r, lin[t].rz, lin[t].rθ),
			 rz = RZLin(model, lin[t].rz),
			 rθ = RθLin(model, lin[t].rθ),
			 v_pr = view(zeros(1,1), 1,1),
			 v_du = view(zeros(1,1), 1,1),
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

function update2!(im_traj::ImplicitTraj, ref_traj::ContactTraj,
	model::ContactDynamicsModel; κ = ref_traj.κ[1],
	)

	H = ref_traj.H

	nq = model.dim.q
	nu = model.dim.u
	nc = model.dim.c
	nb = model.dim.b
	nd = nq + nc + nb
	nz = num_var(model)
	nθ = num_data(model)

	for t = 1:H
		fill!(im_traj.ip[t].κ, κ)
		update!(im_traj.lin[t], model, ref_traj.z[t], ref_traj.θ[t])

		z0  = im_traj.lin[t].z
		θ0  = im_traj.lin[t].θ
		r0  = im_traj.lin[t].r
		rz0 = im_traj.lin[t].rz
		rθ0 = im_traj.lin[t].rθ
		# im_traj.ip[t].r = RLin(model, im_traj.lin[t].z, im_traj.lin[t].θ, im_traj.lin[t].r, im_traj.lin[t].rz, im_traj.lin[t].rθ)#TODO coment
		# im_traj.ip[t].r̄ = RLin(model, im_traj.lin[t].z, im_traj.lin[t].θ, im_traj.lin[t].r, im_traj.lin[t].rz, im_traj.lin[t].rθ)#TODO coment
		# im_traj.ip[t].rz = RZLin(model, im_traj.lin[t].rz)
		# im_traj.ip[t].rθ = RθLin(model, im_traj.lin[t].rθ)
		update!(im_traj.ip[t].r, z0, θ0, r0, rz0, rθ0)# TODO uncommment
		update!(im_traj.ip[t].r̄, z0, θ0, r0, rz0, rθ0)# TODO uncommment
		update!(im_traj.ip[t].rz, rz0)
		update!(im_traj.ip[t].rθ, rθ0)
	end
	return nothing
end

"""
	implicit_dynamics!(im_traj::ImplicitTraj, model::ContactDynamicsModel,
		traj::ContactTraj; κ = traj.κ)
Compute the evaluations and Jacobians of the implicit dynamics on the trajectory 'traj'. The computation is
linearizedimated since it relies on a linearization about a reference trajectory.
"""
function implicit_dynamics!(im_traj::ImplicitTraj, model::ContactDynamicsModel,
	traj; κ = traj.κ) where {T, D}

	# Threads.@threads for t = 1:traj.H
	for t = 1:traj.H
		#@assert abs.(im_traj.lin[t].κ - κ[1]) / κ[1] < 1.0e-5 # check that the κ are consistent between the optimized trajectory (traj)
		# and the linearization (impl.lin).

		# initialized solver
		z_initialize!(im_traj.ip[t].z, model, copy(traj.q[t+2])) #TODO: try alt. schemes
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
