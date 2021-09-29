const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "quadruped_simple", "visuals.jl"))
vis = Visualizer()
open(vis)

s = get_simulation("quadruped_simple", "flat_3D_lc", "flat")
model = s.model
env = s.env
ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/quadruped_simple/gaits/gait2_mit_2.5percent.jld2"),
    load_type = :split_traj_alt))

function check_traj(model::ContactModel, traj::ContactTraj)
    nq = model.dim.q
    nz = num_var(model, env)
    d = [zeros(nq) for t=1:traj.H]
    r = [zeros(nz) for t=1:traj.H]
    for t = 1:ref_traj.H
        k = kinematics(model, traj.q[t+2])
        λ = contact_forces(model, env, traj.γ[t], traj.b[t], traj.q[t+2], k)
		Λ = transpose(J_func(model, env, traj.q[t+2])) * λ
        d[t] = dynamics(model, traj.h, traj.q[t], traj.q[t+1], traj.u[t], traj.w[t], Λ, traj.q[t+2])
        r[t] = residual(model, env, traj.z[t], traj.θ[t], [1e-8])
    end
    return d, r
end
ref_traj
nz = num_var(model, env)
dvio, rvio = check_traj(model, ref_traj)
plot(hcat([v[1:end] for v in dvio]...)')
plot(hcat([v[1:1] for v in dvio]...)')
plot(hcat([v[1:18] for v in rvio]...)')
plot(hcat([v[1:66] for v in rvio]...)')
plot(hcat([v[1:3] for v in rvio]...)')
plot(hcat([v[4:6] for v in rvio]...)')
plot(hcat([v[7:end] for v in rvio]...)')

plot(hcat([v[[1,4,7,10]] for v in ref_traj.u]...)')
plot(hcat([v[[2,5,8,11]] for v in ref_traj.u]...)')
plot(hcat([v[[3,6,9,12]] for v in ref_traj.u]...)')


# time
H = ref_traj.H
h = ref_traj.h
N_sample = 2
H_mpc = 20
h_sim = h / N_sample
H_sim = 300 #4000 #3000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
	q = [Diagonal(1e-1 * [1e-1; 1e0; 1e0; 1e0*ones(3); 3.0 * ones(model.dim.q-6)]) for t = 1:H_mpc],
	# q = [Diagonal(1e-1 * [1e-1; 1e0; 1e0; 1e0*ones(3); [10,3,3]; [10,3,3]; [10,3,3]; [10,3,3]]) for t = 1:H_mpc],
	v = [Diagonal(1e-2 * [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 100.0 * ones(model.dim.q-6)]) for t = 1:H_mpc],
	# v = [Diagonal(1e-2 * [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; [30,3,30]; [30,3,30]; [30,3,30]; [30,3,30]]) for t = 1:H_mpc],
    u = [Diagonal(1e-2 * [1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

get_stride(model, ref_traj)

# nz = num_var(model)
# nθ = num_data(model)
#
# rz0 = zeros(nz,nz)
# rθ0 = zeros(nz,nθ)
# model.res.rz!(rz0, ref_traj.z[1], ref_traj.θ[1])
# model.res.rθ!(rθ0, ref_traj.z[1], ref_traj.θ[1])
# sum(1e10*abs.(rz0) .> 1e5)
# rz0
#
#
# ix, iy1, iy2 = linearization_var_index(model)
# idyn, irst, ibil, ialt = linearization_term_index(model)
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [ix;  ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [iy1; ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [iy2; ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[irst;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[ibil;], [ix; iy1; iy2]])))
#
#
# plot(Gray.(1e10*abs.(rθ0[[idyn; irst; ibil], :])))
# plot(Gray.(1e10*abs.(rθ0[idyn,:])))
# plot(Gray.(1e10*abs.(rθ0[irst,:])))
# plot(Gray.(1e10*abs.(rθ0[ibil,:])))

# # get model
# model_hopper = get_model("hopper_2D")
#
# # get trajectory
# ref_traj_hopper = get_trajectory("hopper_2D", "gait_in_place", load_type = :joint_traj)
#
# nz = num_var(model_hopper)
# nθ = num_data(model_hopper)
#
# rz0 = zeros(nz,nz)
# rθ0 = zeros(nz,nθ)
# model_hopper.res.rz!(rz0, ref_traj_hopper.z[1], ref_traj_hopper.θ[1])
# model_hopper.res.rθ!(rθ0, ref_traj_hopper.z[1], ref_traj_hopper.θ[1])
#
# ix, iy1, iy2 = linearization_var_index(model_hopper)
# idyn, irst, ibil, ialt = linearization_term_index(model_hopper)
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [ix;  ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [iy1; ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [iy2; ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[irst;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[ibil;], [ix; iy1; iy2]])))
#
#
# plot(Gray.(1e10*abs.(rθ0[[idyn; irst; ibil], :])))
# plot(Gray.(1e10*abs.(rθ0[idyn,:])))
# plot(Gray.(1e10*abs.(rθ0[irst,:])))
# plot(Gray.(1e10*abs.(rθ0[ibil,:])))

# function rz!_(model, z, θ, k)
#     f(z) = residual(model, z, θ, [k])
#     return ForwardDiff.jacobian(f, z)
# end
#
#
# rz0 = rz!_(model, ref_traj.z[1], ref_traj.θ[1], 1e-4)
# sum(1e10*abs.(rz0) .> 1e0)
#
#
# ix, iy1, iy2 = linearization_var_index(model)
# idyn, irst, ibil, ialt = linearization_term_index(model)
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [ix;  ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [iy1; ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn; irst; ibil;], [iy2; ]])))
# plot(Gray.(1e10*abs.(rz0[[idyn;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[irst;], [ix; iy1; iy2]])))
# plot(Gray.(1e10*abs.(rz0[[ibil;], [ix; iy1; iy2]])))
#
#
# plot(Gray.(1e10*abs.(rθ0[[idyn; irst; ibil], :])))
# plot(Gray.(1e10*abs.(rθ0[idyn,:])))
# plot(Gray.(1e10*abs.(rθ0[irst,:])))
# plot(Gray.(1e10*abs.(rθ0[ibil,:])))

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        # solver = :ldl_solver,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
		# live_plotting=true
		))

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	# mode = :configurationforce,
	mode = :configuration,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5,
		# max_time = ref_traj.h, # HARD REAL TIME
		),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        altitude_impact_threshold = 0.02,
        altitude_verbose = true,
        ),
	ip_opts = InteriorPointOptions(
		max_iter = 100,
		# verbose = true,
		r_tol = 1.0e-4,
		κ_tol = 1.0e-4,
		diff_sol = true,
		# κ_reg = 1e-3,
		# γ_reg = 1e-1,
		solver = :empty_solver,
		),
    )

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8


sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true),
    )
#
# rz0 = sim.p.im_traj.ip[1].rz
# Ai0 = Array(rz0.S.Ai)
# plot(Gray.(1e10*abs.(Ai0)))


status = ContactControl.simulate!(sim, verbose = true)
# @profiler status = ContactControl.simulate!(sim, verbose = true)

################################################################################
# Timing result
################################################################################
process!(sim)
# Time budget
ref_traj.h
# Time used on average
sim.stats.μ_dt
# Speed ratio
H_sim * h_sim / sum(sim.stats.dt)


# plot_lines!(vis, model, sim.traj.q[1:25:end])
plot_surface!(vis, env, ylims=[0.15, -0.15], xlims=[-1.0, 6.0])
anim = visualize_robot!(vis, model, ref_traj, name=:Ref, sample=1, α=0.5)
anim = visualize_robot!(vis, model, sim.traj, anim=anim, sample=N_sample)
# anim = visualize_robot!(vis, model, sim.traj, anim=anim)
# anim = visualize_force!(vis, model, sim.traj, anim=anim, h=h_sim)



rep_ref_traj = repeat_ref_traj(ref_traj, 10; idx_shift = [1,7,10,13,16])
anim = visualize_robot!(vis, model, rep_ref_traj, name=:Ref, sample=1, α=0.5)



# Display ghosts
t_ghosts = [1, 1333, 2666]
mvis_ghosts = []
for (i,t) in enumerate(t_ghosts)
    α = i/(length(t_ghosts)+1)
    name = Symbol("ghost$i")
    mvis = build_meshrobot!(vis, model, name=name, α=α)
    push!(mvis_ghosts, mvis)
end

for (i,t) in enumerate(t_ghosts)
    name = Symbol("ghost$i")
    set_meshrobot!(vis, mvis_ghosts[i], model, sim.traj.q[t], name=name)
end


# filename = "tablebot_demo"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
