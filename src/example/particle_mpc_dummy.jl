model = get_model("particle", surf = "sinusoidal")
ref_traj = ContactControl.get_trajectory("particle", "sinusoidal2",
    model_name = "particle_sinusoidal")
T = Float64
κ = 1.0e-4
ref_traj.κ[1] = κ
rep_traj = repeat_ref_traj(ref_traj, model, 1, idx_shift = (1:1))
rep_traj_copy = repeat_ref_traj(ref_traj, model, 1, idx_shift = (1:1))

H = rep_traj.H
h = 0.01

N_sample = 3
ref_traj0 = deepcopy(rep_traj)
n_opts0 = NewtonOptions(r_tol=3e-4, κ_init=κ, κ_tol=2κ, solver_inner_iter=5)
m_opts0 = MPCOptions{T}(
            N_sample=N_sample,
            M=500,
            H_mpc=10,
            κ=κ,
            κ_sim=1e-8,
            r_tol_sim=1e-8,
            open_loop_mpc=false,
            w_amp=[-0.0, 0.0, 0.0],
            live_plotting=false)
cost0 = CostFunction(H, model.dim,
    q = [Diagonal(1.0 * ones(model.dim.q))    for t = 1:H],
    u = [Diagonal(1.0e-1 * ones(model.dim.u)) for t = 1:H],
    γ = [Diagonal(1.0e-6 * ones(model.dim.c)) for t = 1:H],
    b = [Diagonal(1.0e-6 * ones(model.dim.b)) for t = 1:H])
core0 = Newton(m_opts0.H_mpc, h, model, cost=cost0, opts=n_opts0)
mpc0 = MPC(model, ref_traj0, m_opts=m_opts0)
@time dummy_mpc(model, core0, mpc0)

qq = []
for q in rep_traj_copy.q
    for i = 1:N_sample
        push!(qq, q)
    end
end

plot(hcat(qq...)[1:3, 1:300]', label = ["x" "y" "z"], color = :black, width = 3.0)
plot!(hcat(Array.(mpc0.q_sim)...)[1:3, 1:300]',label = ["x" "y" "z"], color = :cyan, width = 1.0, legend = :topleft)
