"""
	ImplicitTraj{T,nd}
This structure holds the trajectory of evaluations and Jacobians of the implicit dynamics.
These evaluations and Jacobians are computed using a linear approximation computed around `lin`.

"""
mutable struct ImplicitTraj{T,nd}
	H::Int                               # horizon length
	lin::Vector{LinStep{T}}              # linearization point length=H
	d::Vector{SizedVector{nd,T}}         # dynamics violation  length=H
	δz::Vector{SparseMatrixCSC{T,Int}}   # solution gradient   length=H
end

"""
	Create dummy ImplicitTraj.
"""
function ImplicitTraj(H::Int, model::ContactDynamicsModel)
	nq = model.dim.q
	nc = model.dim.c
	nb = model.dim.b
	nd = nq+nc+nb
	nz = num_var(model)
	nθ = num_data(model)
	lin = [LinStep(model) for k = 1:H]
	d = [zeros(SizedVector{nd}) for k=1:H]
	δz = [spzeros(nz,nθ) for k=1:H]
	return ImplicitTraj{eltype(d[1]),nd}(H,lin,d,δz)
end

"""
	linearization!(model::ContactDynamicsModel, ref_traj::ContactTraj14{T,nq,nu,nw,nc,nb,nz,nθ},
		impl::ImplicitTraj{T,nd}, κ::T=ref_traj.κ) where {T,nq,nu,nw,nc,nb,nz,nθ,nd}
Linearization of the model around the reference trajectory `ref_traj`, the resulting linearization is stored
in `impl.lin`.
"""
function linearization!(model::ContactDynamicsModel, ref_traj::ContactTraj14{T,nq,nu,nw,nc,nb,nz,nθ},
	impl::ImplicitTraj{T,nd}, κ::T=ref_traj.κ) where {T,nq,nu,nw,nc,nb,nz,nθ,nd}
	@assert impl.H == ref_traj.H
	H = ref_traj.H
	for k = 1:H
		# Compute linearization
		# we need to figure out, the correct way to handle ref_traj that are optimized at κ=1e-10,
		# and linearization at κ=1e-3
		impl.lin[k] = LinStep(model, ref_traj.z[k], ref_traj.θ[k], κ)
	end
	return nothing
end

"""
		implicit_dynamics!(model::ContactDynamicsModel, traj::ContactTraj14{T,nq,nu,nw,nc,nb},
	impl::ImplicitTraj{T,nd}; κ::T=traj.κ) where {T,nq,nu,nw,nc,nb,nd}
Compute the evaluations and Jacobians of the implicit dynamics on the trajectory 'traj'. The computation is
approximated since it relies on a linearization about a reference trajectory.
"""
function implicit_dynamics!(model::ContactDynamicsModel, traj::ContactTraj14{T,nq,nu,nw,nc,nb},
	impl::ImplicitTraj{T,nd}; κ::T=traj.κ) where {T,nq,nu,nw,nc,nb,nd}
	@assert impl.H == traj.H
	H = traj.H
	# Compute the implicit dynamics
	# constraint violation (solve linearized contact problem)
	# constraint gradient (implicit function theorem)

	# INITIALIZATION
	ip = interior_point(
		num_var(model),
		num_data(model),
		inequality_indices(model),
		rz = model.spa.rz_sp,
		rθ = model.spa.rθ_sp)
	ip_opts = InteriorPointOptions(
		κ_init=κ,
		κ_tol=κ,
		diff_sol=true)

	for k = 1:H
		@assert abs.(impl.lin[k].κ0 - κ)/κ < 1e-5 # check that the κ are consistent between the optimized trajectory (traj)
		# and the linearization (impl.lin).

		# Define local residual functions
		# residual
		r!(r, z, θ, κ) = r_approx!(impl.lin[k], r, z, θ, κ)
		# residual Jacobian wrt z
		rz!(rz, z, θ, κ) = rz_approx!(impl.lin[k], rz, z, θ, κ)
		# residual Jacobian wrt θ
		rθ!(rθ, z, θ, κ) = rθ_approx!(impl.lin[k], rθ, z, θ, κ)
		# Set the residual functions
		ip.methods.r! = r!
		ip.methods.rz! = rz!
		ip.methods.rθ! = rθ!

		z_initialize!(ip.z, model, traj.q[k+2]) # initialize with our best guess.
		status = interior_point!(ip, ip.z, traj.θ[k]; opts = ip_opts)

		q2 = traj.q[k+2]
		γ1 = traj.γ[k]
		b1 = traj.b[k]

		impl.d[k] = ip.z[1:nq+nc+nb] - [q2; γ1; b1]
		impl.δz[k] = copy(ip.δz)
	end
	return nothing
end
