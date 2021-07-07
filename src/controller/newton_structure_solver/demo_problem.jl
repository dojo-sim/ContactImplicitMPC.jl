# Flamingo
include(joinpath(pwd(), "src/controller/lci_dynamics.jl"))
include(joinpath(pwd(), "src/controller/newton_structure_solver/methods.jl"))

sim = get_simulation("flamingo", "flat_2D_lc", "flat")
model = sim.model
env = sim.env

ref_traj = deepcopy(ContactControl.get_trajectory(model, flat_2D_lc,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

H_mpc = 10
obj_mpc = quadratic_objective(model, H_mpc)#,
    # q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc+2],
    # v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
    # u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(model.dim.u-6); 2; 2]) for t = 1:H_mpc-1])

s = newton_structure_solver(model.dim.q, model.dim.u, H_mpc)
update_objective!(s, obj_mpc)

for t = 1:s.H
	s.qa[t] = ref_traj.q[t] #+ 0.01 * randn(model.dim.q)
	s.qb[t] = ref_traj.q[t+1] #+ 0.01 * randn(model.dim.q)
	t == s.H && continue
	s.u[t] = ref_traj.u[t] + 0.0 * randn(model.dim.u)
	s.ν1[t] = zeros(s.nq)
	s.ν2[t] = zeros(s.nq)
end

for t = 1:s.H+1
	s.q_ref[t] = ref_traj.q[t]
	t >= s.H && continue
	s.u_ref[t] = ref_traj.u[t]
end

ip_opts = eval(interior_point_options(:interior_point))(
			κ_init = 1.0e-4,
			κ_tol = 2.0 * 1.0e-4,
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver)

lci_traj = LCIDynamicsTrajectory(ref_traj, sim,
	ip_type = :interior_point,
	κ = 1.0e-4,
	μ = 0.1,
	opts=ip_opts)

linear_contact_implicit_dynamics!(lci_traj, s.u, s.qa, s.qb, s.H-1)
s.H
s.Aa
s.Ab
s.Ba
s.u
s.rdyn1
s.rdyn2
s.u
s.qa
s.qb
lci_traj
compute_residual!(s, s.u, s.qa, s.qb, s.ν1, s.ν2, lci_traj)


r1 = vcat([[s.rlagu[t]; s.rlagqa[t]; s.rlagqb[t]] for t = 1:s.H-1]...)
r2 = vcat([[s.rdyn1[t]; s.rdyn2[t]] for t = 1:s.H-1]...)
norm(r1,1)
norm(r2,1)
res1 = [r1; r2]
norm(res1, 1)
residual_norm(s)

factorize!(s)
ContactControl.solve!(s)

α = 0.25
step!(s, α)

compute_residual!(s, s.u_cand, s.qa_cand, s.qb_cand, s.ν1_cand, s.ν2_cand, lci_traj)

s1 = vcat([[s.rlagu[t]; s.rlagqa[t]; s.rlagqb[t]] for t = 1:s.H-1]...)
s2 = vcat([[s.rdyn1[t]; s.rdyn2[t]] for t = 1:s.H-1]...)
norm(s1,1)
norm(s2,1)
res2 = [s1; s2]
norm(res2, 1)
residual_norm(s)


accept_step!(s)
