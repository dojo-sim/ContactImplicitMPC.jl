# Newton solver options
@with_kw mutable struct NewtonOptions{T}
    r_tol::T = 1.0e-5            # primal dual residual tolerance
	max_iter::Int = 10           # maximum number of primal dual iter
    max_time::T = 10000.0        # maximum time spent in the Newton solver
    β_init::T = 1.0e-5           # initial dual regularization
    live_plotting::Bool = false  # visualize the trajectory during the solve
    threads::Bool=false
    verbose::Bool = false
    solver::Symbol = :lu_solver  # lu_sparse_solver
end

abstract type NewtonIndices end
abstract type NewtonResidual end
abstract type NewtonJacobian end

mutable struct Newton{T,nq,nu,nw,nc,nb,nz,nθ,nν,NJ,NR,NI,O,LS,NV}
    jac::NJ                         # NewtonJacobian
    res::NR                        # residual
    res_cand::NR                     # candidate residual
    Δ::NR                          # step direction in the Newton solve, it contains: q2-qH+1, u1-uH, γ1-γH, b1-bH, λd1-λdH
    ν::Vector{NV}          # implicit dynamics lagrange multiplier
    ν_cand::Vector{NV}         # candidate implicit dynamics lagrange multiplier
    traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}                           # optimized trajectory
    traj_cand::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ}      # trial trajectory used in line search
    Δq::Vector{Vector{T}}#::Vector{SizedArray{Tuple{nq},T,1,1}}         # difference between the traj and ref_traj
    Δu::Vector{Vector{T}}         # difference between the traj and ref_traj
    Δγ::Vector{Vector{T}}         # difference between the traj and ref_traj
    Δb::Vector{Vector{T}}         # difference between the traj and ref_traj
    ind::NI                                 # indices of a one-time-step block
    obj::O                              # obj function
    solver::LS
    β::T
    opts::NewtonOptions{T}                          # Newton solver options
end

function Newton(s::Simulation{T}, H::Int, h::T,
    traj::ContactTraj{T}, im_traj::ImplicitTrajectory{T};
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
    window = collect(1:(H + 2))
    implicit_dynamics!(im_traj, traj, window=window, threads=opts.threads) #@@@
    jacobian!(jac, im_traj, obj, H, opts.β_init, window)

    res = NewtonResidual(model, env, H, mode = mode)
    res_cand = NewtonResidual(model, env, H, mode = mode)

    Δ = NewtonResidual(model, env, H, mode = mode)

    ν = [zeros(SizedVector{ind.nd,T}) for t = 1:H]
    ν_cand = deepcopy(ν)

    traj = contact_trajectory(model, env, H, h, κ=κ)
    traj_cand = contact_trajectory(model, env, H, h, κ=κ)

    Δq  = [zeros(nq) for t = 1:H]

    Δu  = [zeros(nu) for t = 1:H]
    Δγ  = [zeros(nc) for t = 1:H]
    Δb  = [zeros(nb) for t = 1:H]

    # regularization
    β = copy(opts.β_init)

    # linear solver
    solver = eval(opts.solver)(jac.R)

    return Newton{T,nq,nu,nw,nc,nb,nz,nθ,nd,typeof.([jac,res,ind,obj,solver,ν[1]])...}(
        jac, res, res_cand, Δ, ν, ν_cand, traj, traj_cand,
        Δq, Δu, Δγ, Δb, ind, obj, solver, β, opts)
end

function delta!(Δx::SizedArray{Tuple{nx},T,1,1}, x, x_ref) where {nx,T}
    Δx .= x
    Δx .-= x_ref
    return nothing
end

function delta!(Δx::Vector{T}, x, x_ref) where T
    Δx .= x
    Δx .-= x_ref
    return nothing
end

function copy_traj!(traj::ContactTraj, traj_cand::ContactTraj, H::Int)
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

function reset!(core::Newton, ref_traj::ContactTraj,
    q0::Vector{T}, q1::Vector{T};
    window=collect(1:(core.traj.H + 2)),
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

        # TODO: not sure this is correct
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
    window::Vector{Int},
    im_traj::ImplicitTrajectory{T},
    ref_traj::ContactTraj{T};
    warm_start::Bool=false) where T

    elapsed_time = 0.0

    elapsed_time += @elapsed begin
        # reset solver
        reset!(core, ref_traj, q0, q1,
            window=window,
            warm_start=warm_start)
    end
    elapsed_time >= core.opts.max_time && (return nothing)

    elapsed_time += @elapsed begin
        # Compute implicit dynamics about traj
        implicit_dynamics!(im_traj, core.traj, window=window, threads=core.opts.threads)
    end
    elapsed_time >= core.opts.max_time && (return nothing)

    elapsed_time += @elapsed begin
        # Compute residual
        residual!(core.res, core, core.ν, im_traj, core.traj, ref_traj, window)
        r_norm = norm(core.res.r, 1)
    end
    elapsed_time >= core.opts.max_time && (return nothing)

    for l = 1:core.opts.max_iter
		# elapsed_time >= core.opts.max_time && break
		# elapsed_time += @elapsed begin
        # check convergence
        r_norm / length(core.res.r) < core.opts.r_tol && break

        elapsed_time += @elapsed begin
            # Compute NewtonJacobian
            jacobian!(core.jac, im_traj, core.obj, core.traj.H, core.β, window, update_hessian=true)
        end
        elapsed_time >= core.opts.max_time && (return nothing)
        # jacobian!(core.jac, im_traj, core.obj, core.traj.H, core.β, update_hessian=false)
        # @show diag(core.jac.R[1:6,1:6])

        elapsed_time += @elapsed begin
            # Compute Search Direction
            linear_solve!(core.solver, core.Δ.r, core.jac.R, core.res.r)
        end
        elapsed_time >= core.opts.max_time && (return nothing)

        # line search the step direction
        α = 1.0
        iter = 0

        elapsed_time += @elapsed begin
            # candidate step
            update_traj!(core.traj_cand, core.traj, core.ν_cand, core.ν, core.Δ, α)
        end
        elapsed_time >= core.opts.max_time && (return nothing)

        elapsed_time += @elapsed begin
            # Compute implicit dynamics for candidate
            implicit_dynamics!(im_traj, core.traj_cand, window=window, threads=core.opts.threads)
        end
        elapsed_time >= core.opts.max_time && (return nothing)

        elapsed_time += @elapsed begin
            # Compute residual for candidate
            residual!(core.res_cand, core, core.ν_cand, im_traj, core.traj_cand, ref_traj, window)
            r_cand_norm = norm(core.res_cand.r, 1)
        end
        elapsed_time >= core.opts.max_time && (return nothing)

        while r_cand_norm^2.0 >= (1.0 - 0.001 * α) * r_norm^2.0
            α = 0.5 * α

            iter += 1
            if iter > 6
                break
            end

            elapsed_time += @elapsed begin
                update_traj!(core.traj_cand, core.traj, core.ν_cand, core.ν, core.Δ, α)
            end
            elapsed_time >= core.opts.max_time && (return nothing)

            elapsed_time += @elapsed begin
                # Compute implicit dynamics about trial_traj
                implicit_dynamics!(im_traj, core.traj_cand, window=window, threads=core.opts.threads)
            end
            elapsed_time >= core.opts.max_time && (return nothing)

            elapsed_time += @elapsed begin
                residual!(core.res_cand, core, core.ν_cand, im_traj, core.traj_cand, ref_traj, window)
                r_cand_norm = norm(core.res_cand.r, 1)
            end
            elapsed_time >= core.opts.max_time && (return nothing)
        end

        elapsed_time += @elapsed begin
            # update
            update_traj!(core.traj, core.traj, core.ν, core.ν, core.Δ, α)
            core.res.r .= core.res_cand.r
            r_norm = r_cand_norm
        end
        elapsed_time >= core.opts.max_time && (return nothing)

        # regularization update
        iter > 6 ? (core.β = min(core.β * 1.3, 1.0e2)) : (core.β = max(1.0e1, core.β / 1.3))

        # print
        core.opts.verbose && print_status(core, im_traj, elapsed_time, α, l)
		# end
    end

    return nothing
end

function print_status(core::Newton, im_traj, elapsed_time, α, l)
     # print status
     println(
     "     t: ", scn(elapsed_time, digits=0),
     "     r̄: ", scn(norm(core.res_cand.r, 1) / length(core.res_cand.r), digits=0),
     "     r: ", scn(norm(core.res.r, 1) / length(core.res.r), digits=0),
     "     Δ: ", scn(norm(core.Δ.r, 1) / length(core.Δ.r), digits=0),
     "     α: ", -Int(round(log(α))),
     "     l: ", l,
     "     κ: ", scn(im_traj.ip[1].κ[1], digits = 0)
     )
end

