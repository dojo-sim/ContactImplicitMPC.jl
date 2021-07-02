sim = get_simulation("flamingo", "flat_2D_lc", "flat")
model = sim.model
env = sim.env

H_mpc = 10
s = newton_structure_solver(model.dim.q, model.dim.u, H_mpc)
obj_mpc = quadratic_objective(model, H_mpc)

update_objective!(s, obj_mpc)

ref_traj = deepcopy(ContactControl.get_trajectory(model, flat_2D_lc,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

traj = deepcopy(ref_traj)

ip_opts = eval(interior_point_options(:interior_point))(
			κ_init = 1.0e-4,
			κ_tol = 2.0 * 1.0e-4,
			r_tol = 1.0e-8,
			diff_sol = true,
			solver = :empty_solver)
im_traj = ImplicitTraj(traj, sim,
	ip_type = :interior_point,
	κ = 1.0e-4,
	opts=ip_opts,
	mode = :configuration)

implicit_dynamics!(im_traj, sim, traj, κ = [1.0e-4])

im_traj
