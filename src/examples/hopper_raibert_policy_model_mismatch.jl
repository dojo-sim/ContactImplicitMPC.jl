include(joinpath(@__DIR__, "..", "dynamics", "hopper_2D", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)

# nominal model
include(joinpath(pwd(), "src/dynamics/hopper_2D/model.jl"))
model_nom = get_model("hopper_2D")

model_1 = Hopper2D(Dimensions(nq, nu, nw, nc, nb),
			   1.0 * mb, ml, 1.0 * Jb, Jl,
			   μ_world, μ_joint, g,
			   BaseMethods(), DynamicsMethods(), ContactMethods(),
			   ResidualMethods(), ResidualMethods(),
			   SparseStructure(spzeros(0, 0), spzeros(0, 0)),
			   SVector{4}(zeros(4)),
			   environment_2D_flat())

expr_base = generate_base_expressions_analytical(model_1)
instantiate_base!(model_1.base, expr_base)
expr_dyn = generate_dynamics_expressions(model_1)
instantiate_dynamics!(model_1.dyn, expr_dyn)
expr_res, rz_sp, rθ_sp = generate_residual_expressions(model_1)
instantiate_residual!(model_1.res, expr_res)
model_1.spa.rz_sp = rz_sp
model_1.spa.rθ_sp = rθ_sp

nq = model_nom.dim.q
nu = model_nom.dim.u
nc = model_nom.dim.c
nb = model_nom.dim.b

# time
H = 92#traj.H
h = 0.01#traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 5000

# Select the initial speed
v0 = 0.2
Tstance = 0.13 # measure using hop-in-place gait
Tflight = 0.62 # measure using hop-in-place gait
include(joinpath(pwd(), "src/controller/raibert_policy.jl"))
p = raibert_policy(model_nom, v0=v0, Tstance=Tstance, Tflight=Tflight, h=h)

off0 = SVector{nq,T}([0.0, 0.5, 0.0, 0.0])
off1 = SVector{nq,T}([v0*h_sim, 0.5, 0.0, 0.0])
q_ref = SVector{nq,T}([0.0, 0.5, 0.0, 0.5])
q0_sim = copy(q_ref) + off0
q1_sim = copy(q_ref) + off1

sim = simulator(model_1, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
	uL = [-1.0, -0.3], uU = [1.0, 0.3],
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true)
    )

@time status = simulate!(sim)

anim = visualize_robot!(vis, sim.model, sim.traj, sample=20)

plot(hcat(sim.traj.u...)[1:2, 1:10:end]')
