# # MPC options
# @with_kw mutable struct MPCOptions{T}
# 	κ::T = 1.0e-4
# 	ip_max_time::T = 60.0     # maximum time allowed for an InteriorPoint solve
#     live_plotting::Bool = false # Use the live plotting tool to debug
# 	altitude::Bool = false
# end
#
# mutable struct MPC{T,nq}
#     ref_traj::ContactTraj
#     im_traj::ImplicitTraj
#     stride::SizedArray{Tuple{nq},T,1,1,Vector{T}} # change in q required to loop the ref trajectories
# 	altitude::Vector{T} # TODO: static arrays
# 	opts::MPCOptions{T}
# end
#
# function MPC(model::ContactDynamicsModel, ref_traj::ContactTraj{T,nq,nu,nw,nc,nb,nz,nθ};
#         opts::MPCOptions = MPCOptions()) where {T,nq,nu,nw,nc,nb,nz,nθ}
#
#     im_traj = ImplicitTraj(ref_traj, model, κ = opts.κ, max_time = opts.ip_max_time)
#     stride = get_stride(model, ref_traj)
# 	altitude = zeros(model.dim.c)
#
#     return MPC{T,nq}(ref_traj, im_traj, stride, altitude, opts)
# end

function rot_n_stride!(traj::ContactTraj, stride::SizedArray)
    rotate!(traj)
    mpc_stride!(traj, stride)
    return nothing
end

function rotate!(traj::ContactTraj)
    H = traj.H
	# TODO: remove copy
    q1 = copy(traj.q[1])
    u1 = copy(traj.u[1])
    w1 = copy(traj.w[1])
    γ1 = copy(traj.γ[1])
    b1 = copy(traj.b[1])
    z1 = copy(traj.z[1])
    θ1 = copy(traj.θ[1])

    for t = 1:H+1
        traj.q[t] .= traj.q[t+1]
    end

    for t = 1:H-1
        traj.u[t] .= traj.u[t+1]
        traj.w[t] .= traj.w[t+1]
        traj.γ[t] .= traj.γ[t+1]
        traj.b[t] .= traj.b[t+1]
        traj.z[t] .= traj.z[t+1]
        traj.θ[t] .= traj.θ[t+1]
    end

    traj.q[end] .= q1
    traj.u[end] .= u1
    traj.w[end] .= w1
    traj.γ[end] .= γ1
    traj.b[end] .= b1
    traj.z[end] .= z1
    traj.θ[end] .= θ1

    return nothing
end

"""
    Update the last two cofiguaraton to be equal to the fisrt two up to an constant offset.
"""
function mpc_stride!(traj::ContactTraj, stride::SizedArray)
    H = traj.H

    for t = H+1:H+2
        traj.q[t] .= deepcopy(traj.q[t-H] + stride)
        update_z!(traj, t - 2)
        update_θ!(traj, t - 2)
        update_θ!(traj, t - 3)
    end

    return nothing
end

function get_stride(model::ContactDynamicsModel, traj::ContactTraj) #TODO: dispatch over environment / model
    stride = zeros(SizedVector{model.dim.q})
    stride[1] = traj.q[end-1][1] - traj.q[1][1]
    return stride
end

function update_altitude!(alt, q)
	return nothing
end

# function dummy_mpc(model::ContactDynamicsModel, core::Newton, mpc::MPC; verbose::Bool=false)
#     elap = 0.0
#     # Unpack
#     # impl = mpc.impl
#     ref_traj = mpc.ref_traj
#     m_opts = mpc.m_opts
#     cost = core.cost
#     opts = core.opts
#     N_sample = mpc.m_opts.N_sample
#     h_sim = mpc.ref_traj.h/N_sample
#     nq = model.dim.q
#
#     # Initialization
#     q0 = deepcopy(ref_traj.q[1])
#     q1 = deepcopy(ref_traj.q[2])
#     q0_sim = SVector{nq}(deepcopy(q1 + (q0 - q1)/N_sample))
#     q1_sim = SVector{nq}(deepcopy(q1))
#     push!(mpc.q_sim, [q0_sim, q1_sim]...)
#
#     for l = 1:m_opts.M
#
#         # Get control
#         if m_opts.open_loop_mpc
#             u_zoh  = SVector{model.dim.u}.([deepcopy(ref_traj.u[1])/N_sample for j=1:N_sample])
#         else
#             warm_start = l > 1
#             m_opts.altitude && update_altitude!(mpc.alt, mpc.q_sim[end])
#             elap += 0.0 * @elapsed update!(mpc.impl, ref_traj, model, mpc.altitude, κ=m_opts.κ)
#             elap += @elapsed newton_solve!(core, model, mpc.impl, ref_traj; verbose=verbose, warm_start=warm_start, q0=q0, q1=q1)
#             u_zoh  = SVector{model.dim.u}.([deepcopy(core.traj.u[1])/N_sample for j=1:N_sample])
#         end
#         # Get disturbances
#         # Apply control and rollout dynamics
#         sim = simulator(model, SVector{model.dim.q}(deepcopy(mpc.q_sim[end-1])), SVector{model.dim.q}(deepcopy(mpc.q_sim[end])), h_sim, N_sample,
#             p = open_loop_policy(u_zoh),
#             d = random_disturbances(model, m_opts.w_amp/N_sample, N_sample, h_sim),
#             r! = model.res.r!,
#             rz! = model.res.rz!,
#             rθ! = model.res.rθ!,
#             rz = model.spa.rz_sp,
#             rθ = model.spa.rθ_sp,
#             ip_opts = InteriorPointOptions(r_tol = m_opts.r_tol_sim, κ_tol = 2m_opts.κ_sim, κ_init = m_opts.κ_sim),
#             sim_opts = SimulatorOptions(warmstart = true)
#             )
#         simulate!(sim)
#         ########################
#         if m_opts.live_plotting
#             plt = plot(layout=grid(2,1,heights=[0.7, 0.3], figsize=[(1000, 1000),(400,400)]), legend=false, xlims=(0,20))
#             plot!(plt[1,1], hcat(Vector.(core.traj.q)...)', color=:blue, linewidth=1.0)
#             plot!(plt[1,1], hcat(Vector.(ref_traj.q)...)', linestyle=:dot, color=:red, linewidth=3.0)
#
#             scatter!((2-1/N_sample)ones(model.dim.q), sim.traj.q[1], markersize=8.0, color=:lightgreen)
#             scatter!(2*ones(model.dim.q), sim.traj.q[2], markersize=8.0, color=:lightgreen)
#             scatter!(plt[1,1], 3*ones(model.dim.q), sim.traj.q[end], markersize=8.0, color=:lightgreen)
#
#             scatter!(plt[1,1], 1*ones(model.dim.q), q0, markersize=6.0, color=:blue)
#             scatter!(plt[1,1], 2*ones(model.dim.q), q1, markersize=6.0, color=:blue)
#
#             scatter!(plt[1,1], 1*ones(model.dim.q), ref_traj.q[1], markersize=4.0, color=:red)
#             scatter!(plt[1,1], 2*ones(model.dim.q), ref_traj.q[2], markersize=4.0, color=:red)
#             scatter!(plt[1,1], 3*ones(model.dim.q), ref_traj.q[3], markersize=4.0, color=:red)
#
#             plot!(plt[2,1], hcat(Vector.(core.traj.u)...)', color=:blue, linewidth=1.0)
#             plot!(plt[2,1], hcat(Vector.(ref_traj.u)...)', linestyle=:dot, color=:red, linewidth=3.0)
#             display(plt)
#         end
#         ########################
#         push!(mpc.q_sim, deepcopy(sim.traj.q[3:end])...)
#         push!(mpc.u_sim, deepcopy(sim.traj.u)...)
#         push!(mpc.w_sim, deepcopy(sim.traj.w)...)
#         push!(mpc.γ_sim, deepcopy(sim.traj.γ)...)
#         push!(mpc.b_sim, deepcopy(sim.traj.b)...)
#         ########################
#         q0 = deepcopy(q1)
#         q1 = deepcopy(sim.traj.q[end])
#         rot_n_stride!(ref_traj, mpc.q_stride)
#     end
#     @show elap
#     return nothing
# end

function live_plotting(model, ref_traj, sim_traj, newton)
	plt = plot(layout=grid(2,1,heights=[0.7, 0.3], figsize=[(1000, 1000),(400,400)]), legend=false, xlims=(0,20))
	plot!(plt[1,1], hcat(Vector.(neton.traj.q)...)', color=:blue, linewidth=1.0)
	plot!(plt[1,1], hcat(Vector.(ref_traj.q)...)', linestyle=:dot, color=:red, linewidth=3.0)

	scatter!((2-1/N_sample)ones(model.dim.q), sim_traj.q[1], markersize=8.0, color=:lightgreen)
	scatter!(2*ones(model.dim.q), sim_traj.q[2], markersize=8.0, color=:lightgreen)
	scatter!(plt[1,1], 3*ones(model.dim.q), sim_traj.q[end], markersize=8.0, color=:lightgreen)

	# scatter!(plt[1,1], 1*ones(model.dim.q), q0, markersize=6.0, color=:blue)
	# scatter!(plt[1,1], 2*ones(model.dim.q), q1, markersize=6.0, color=:blue)

	scatter!(plt[1,1], 1*ones(model.dim.q), ref_traj.q[1], markersize=4.0, color=:red)
	scatter!(plt[1,1], 2*ones(model.dim.q), ref_traj.q[2], markersize=4.0, color=:red)
	scatter!(plt[1,1], 3*ones(model.dim.q), ref_traj.q[3], markersize=4.0, color=:red)

	plot!(plt[2,1], hcat(Vector.(newton.traj.u)...)', color=:blue, linewidth=1.0)
	plot!(plt[2,1], hcat(Vector.(ref_traj.u)...)', linestyle=:dot, color=:red, linewidth=3.0)
	display(plt)
end
