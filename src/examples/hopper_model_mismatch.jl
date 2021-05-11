include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# simulation model
include(joinpath(pwd(), "src/dynamics/hopper_2D/model.jl"))

model_sim = Hopper2D(Dimensions(nq, nu, nw, nc, nb),
			   1.2 * mb, ml, 1.2 * Jb, Jl,
			   μ_world, μ_joint, g,
			   BaseMethods(), DynamicsMethods(), ContactMethods(),
			   ResidualMethods(), ResidualMethods(),
			   SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			   SVector{4}(zeros(4)),
			   environment_2D_flat())

expr_base = generate_base_expressions_analytical(model_sim)
instantiate_base!(model_sim.base, expr_base)
expr_dyn = generate_dynamics_expressions(model_sim)
instantiate_dynamics!(model_sim.dyn, expr_dyn)
expr_res, rz_sp, rθ_sp = generate_residual_expressions(model_sim)
instantiate_residual!(model_sim.res, expr_res)
model_sim.spa.rz_sp = rz_sp
model_sim.spa.rθ_sp = rθ_sp

# get hopper model
model = get_model("hopper_2D")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
# ref_traj = get_trajectory("hopper_2D", "gait_in_place", load_type=:joint_traj)
ref_traj = get_trajectory("hopper_2D", "gait_forward", load_type=:joint_traj)
ref_traj_copy = deepcopy(ref_traj)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 5000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(H_mpc, model.dim,
    q = [Diagonal(1.0e-1 * [1,3,1,3])   for t = 1:H_mpc],
    u = [Diagonal(1.0e-0 * [1e-2, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
        altitude_update = true,
        altitude_impact_threshold = 0.5,
        altitude_verbose = true,
        )
    )


q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactControl.simulator(model_sim, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	uL = [-0.25, -0.75], uU = [0.25, 0.75],
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))

@time status = ContactControl.simulate!(sim)

(sim.traj.q[end][1] - sim.traj.q[1][1]) / (h_sim * H_sim)
anim = visualize_robot!(vis, sim.model, sim.traj, sample=20)

plot(hcat(sim.traj.u...)[1:2, 1:10:end]')
