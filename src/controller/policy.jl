"""
    linearized model-predictive control policy
"""
mutable struct LinearizedMPC <: Policy
   mpc
   core
   model
   t_prev
   q0
   N_sample
   verbose
   warmstart
   cnt
end

function linearized_mpc_policy(ref_traj, model, cost, H_mpc, M, N_sample;
   verbose = false,
   warmstart = false,
   κ = 1.0e-4,
   n_opts = NewtonOptions(
      r_tol=3e-4,
      κ_init=κ,
      κ_tol=2κ,
      solver_inner_iter=5),
   m_opts = MPCOptions{Float64}(
      N_sample=N_sample,
      M=M,
      H_mpc=H_mpc,
      κ=κ,
      κ_sim=1e-8,
      r_tol_sim=1e-8,
      open_loop_mpc=false,
      w_amp=[0.0, -0.0],
      ip_max_time=100.0,
      live_plotting=false))
   ref_traj0 = deepcopy(ref_traj)

   core0 = Newton(m_opts.H_mpc, ref_traj0.h, model, cost=cost, opts=n_opts)
   mpc0 = MPC(model, ref_traj0, m_opts=m_opts)

   LinearizedMPC(mpc0, core0, model, -ref_traj0.h, copy(ref_traj0.q[1]), N_sample, verbose, warmstart, N_sample)
end


function policy(p::LinearizedMPC, x, traj, t)
    # @show p.cnt
    # @show t
    # @show "$(p.t_prev + p.mpc.ref_traj.h)"
    # @show t >= p.t_prev + p.mpc.ref_traj.h
    if p.cnt == p.N_sample
       # p.mpc.impl = ImplicitTraj(p.mpc.ref_traj, p.model, κ=p.mpc.m_opts.κ, max_time=p.mpc.m_opts.ip_max_time)
       update!(p.mpc.impl, p.mpc.ref_traj, p.model, κ=p.mpc.m_opts.κ, max_time=p.mpc.m_opts.ip_max_time)
       newton_solve!(p.core, p.model, p.mpc.impl, p.mpc.ref_traj; verbose=p.verbose, warm_start= t > 0.0, q0=copy(p.q0), q1=copy(x))
       rot_n_stride!(p.mpc.ref_traj, p.mpc.q_stride)
       p.q0 .= copy(x)
       # p.t_prev = copy(t)
       p.cnt = 0
    end
    p.cnt += 1
    return p.core.traj.u[1] / p.N_sample # rescale output
end
