# Newton solver options
@with_kw mutable struct NewtonOptions{T}
    r_tol::T = 1.0e-5            # primal dual residual tolerance
	max_iter::Int = 10           # maximum number of primal dual iter
    max_time::T = 10000.0        # maximum time spent in the Newton solver
    β_init::T = 1.0e-5           # initial dual regularization
    live_plotting::Bool = false  # visualize the trajectory during the solve
    verbose::Bool = false
    solver::Symbol = :lu_solver
end

abstract type NewtonIndices end
abstract type NewtonResidual end
abstract type NewtonJacobian end

mutable struct Newton{T,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS}
    jac::NJ                         # NewtonJacobian
    res::NR                        # residual
    res_cand::NR                     # candidate residual
    Δ::NR                          # step direction in the Newton solve, it contains: q2-qH+1, u1-uH, γ1-γH, b1-bH, λd1-λdH
    ν::Vector{SizedArray{Tuple{nν},T,1,1}}          # implicit dynamics lagrange multiplier
    ν_cand::Vector{SizedArray{Tuple{nν},T,1,1}}         # candidate implicit dynamics lagrange multiplier
    traj::ContactTrajectory{T,nq,nu,nw,nc,nb,nz,nθ}                           # optimized trajectory
    traj_cand::ContactTrajectory{T,nq,nu,nw,nc,nb,nz,nθ}      # trial trajectory used in line search
    Δq::Vector{SizedArray{Tuple{nq},T,1,1}}         # difference between the traj and ref_traj
    Δu::Vector{SizedArray{Tuple{nu},T,1,1}}         # difference between the traj and ref_traj
    Δγ::Vector{SizedArray{Tuple{nc},T,1,1}}         # difference between the traj and ref_traj
    Δb::Vector{SizedArray{Tuple{nb},T,1,1}}         # difference between the traj and ref_traj
    ind::NI                                 # indices of a one-time-step block
    obj::O                              # obj function
    solver::LS
    β::T
    opts::NewtonOptions{T}                          # Newton solver options
end

function Newton(s::Simulation{T}, H::Int, h::T,
    traj::ContactTrajectory{T}, im_traj::ImplicitTrajectory{T};
    obj::Objective = TrackingObjective(s.model, s.env, H),
    opts::NewtonOptions{T} = NewtonOptions(), κ::T=im_traj.ip[1].κ[1]) where T

    model = s.model
    env = s.env

    mode = im_traj.mode

    ind = NewtonIndices(model, env, H, mode = mode)

    nq = model.nq
    nu = model.nu
    nw = model.nw
    nc = model.nc
    nb = nc * friction_dim(env)
    nz = num_var(model, env)
    nθ = num_data(model)
    nd = ind.nd

    jac = NewtonJacobian(model, env, H, mode = mode)

    # precompute Jacobian for pre-factorization
    implicit_dynamics!(im_traj, s, traj, κ = κ) #@@@

    jacobian!(jac, im_traj, obj, H, opts.β_init)

    res = NewtonResidual(model, env, H, mode = mode)
    res_cand = NewtonResidual(model, env, H, mode = mode)

    Δ = NewtonResidual(model, env, H, mode = mode)

    ν = [zeros(SizedVector{ind.nd,T}) for t = 1:H]
    ν_cand = deepcopy(ν)

    traj = contact_trajectory(model, env, H, h, κ=κ)
    traj_cand = contact_trajectory(model, env, H, h, κ=κ)

    Δq  = [zeros(SizedVector{nq,T}) for t = 1:H]
    Δu  = [zeros(SizedVector{nu,T}) for t = 1:H]
    Δγ  = [zeros(SizedVector{nc,T}) for t = 1:H]
    Δb  = [zeros(SizedVector{nb,T}) for t = 1:H]

    # regularization
    β = copy(opts.β_init)

    # linear solver
    solver = eval(opts.solver)(jac.R)

    return Newton{T,nq,nu,nw,nc,nb,nz,nθ,nd,typeof.([jac,res,ind,obj,solver])...}(
        jac, res, res_cand, Δ, ν, ν_cand, traj, traj_cand,
        Δq, Δu, Δγ, Δb, ind, obj, solver, β, opts)
end

function delta!(Δx::SizedArray{Tuple{nx},T,1,1}, x, x_ref) where {nx,T}
    Δx .= x
    Δx .-= x_ref
    return nothing
end

#TODO: add minus function

function copy_traj!(traj::ContactTrajectory, traj_cand::ContactTrajectory, H::Int)
    Ht = traj.H
    Hs = traj_cand.H # MAYBE BREAKING TEST

    @assert Hs >= H
    @assert Ht >= H

    traj.κ .= traj_cand.κ

    for t in eachindex(1:H + 2)
        traj.q[t] .= traj_cand.q[t]
    end

    for t in eachindex(1:H)
        traj.u[t] .= traj_cand.u[t]
        traj.w[t] .= traj_cand.w[t]
        traj.γ[t] .= traj_cand.γ[t]
        traj.b[t] .= traj_cand.b[t]
        traj.z[t] .= traj_cand.z[t]
        traj.θ[t] .= traj_cand.θ[t]
    end

    return nothing
end

function reset!(core::Newton, ref_traj::ContactTrajectory,
    q0::Vector{T}, q1::Vector{T};
    warm_start::Bool = false) where T

    # H = ref_traj.H
    H_mpc = core.traj.H
    opts = core.opts

    # Reset β value
    core.β = opts.β_init

    if !warm_start
		# Reset duals
		for t = 1:H_mpc
			fill!(core.ν[t], 0.0)
			fill!(core.ν_cand[t], 0.0)
		end

		# Set up trajectory
        copy_traj!(core.traj, ref_traj, core.traj.H)
	end

    core.traj.q[1] .= q0
    core.traj.q[2] .= q1

    update_θ!(core.traj, 1)
    update_θ!(core.traj, 2)

    # initialized residual Jacobian
    initialize_jacobian!(core.jac, core.obj, core.traj.H)

	# Set up traj cand
    copy_traj!(core.traj_cand, core.traj, core.traj.H)

	return nothing
end

function newton_solve!(
    core::Newton{T},
    s::Simulation{T},
    q0::Vector{T},
    q1::Vector{T},
    im_traj::ImplicitTrajectory{T},
    ref_traj::ContactTrajectory{T};
    warm_start::Bool=false) where T

    # reset solver 
    reset!(core, ref_traj, q0, q1, warm_start=warm_start)
    
    # Compute implicit dynamics about traj
	implicit_dynamics!(im_traj, s, core.traj, κ = im_traj.ip[1].κ[1])
    
    # return nothing
    # Compute residual
    residual!(core.res, core, core.ν, im_traj, core.traj, ref_traj)

    r_norm = norm(core.res.r, 1)
	elapsed_time = 0.0

    for l = 1:core.opts.max_iter
		elapsed_time >= core.opts.max_time && break
		# println("iter:", l, "  elapsed_time:", elapsed_time)
		elapsed_time += @elapsed begin
	        # check convergence
	        r_norm / length(core.res.r) < core.opts.r_tol && break
	        # Compute NewtonJacobian
	        update_jacobian!(core.jac, im_traj, core.obj, core.traj.H, core.β)

	        # Compute Search Direction
	        linear_solve!(core.solver, core.Δ.r, core.jac.R, core.res.r)

	        # line search the step direction
	        α = 1.0
	        iter = 0

	        # candidate step
	        update_traj!(core.traj_cand, core.traj, core.ν_cand, core.ν, core.Δ, α)

	        # Compute implicit dynamics for candidate
			implicit_dynamics!(im_traj, s, core.traj_cand, κ = im_traj.ip[1].κ[1])

	        # Compute residual for candidate
	        residual!(core.res_cand, core, core.ν_cand, im_traj, core.traj_cand, ref_traj)
	        r_cand_norm = norm(core.res_cand.r, 1)

	        while r_cand_norm^2.0 >= (1.0 - 0.001 * α) * r_norm^2.0
	            α = 0.5 * α

	            iter += 1
	            if iter > 6
	                break
	            end

	            update_traj!(core.traj_cand, core.traj, core.ν_cand, core.ν, core.Δ, α)

	            # Compute implicit dynamics about trial_traj
				implicit_dynamics!(im_traj, s, core.traj_cand, κ = im_traj.ip[1].κ[1])

	            residual!(core.res_cand, core, core.ν_cand, im_traj, core.traj_cand, ref_traj)
	            r_cand_norm = norm(core.res_cand.r, 1)
	        end

	        # update
	        update_traj!(core.traj, core.traj, core.ν, core.ν, core.Δ, α)
	        core.res.r .= core.res_cand.r
	        r_norm = r_cand_norm
	        # println("lafter = ", l, "  norm = ", r_norm / length(core.res.r)) #@@@

	        # regularization update
	        if iter > 6
	            core.β = min(core.β * 1.3, 1.0e2)
	        else
	            core.β = max(1.0e1, core.β / 1.3)
	        end

	        # print status
	        core.opts.verbose && println(" l: ", l ,
					"     t: ", scn(elapsed_time, digits = 0),
	                "     r̄: ", scn(norm(core.res_cand.r, 1) / length(core.res_cand.r), digits = 0),
	                "     r: ", scn(norm(core.res.r, 1) / length(core.res.r), digits = 0),
	                "     Δ: ", scn(norm(core.Δ.r, 1) / length(core.Δ.r), digits = 0),
	                "     α: ", -Int(round(log(α))),
	                "     κ: ", scn(im_traj.ip[1].κ[1], digits = 0))
		end
    end

    return nothing
end
