const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "quadrupedlinear", "visuals.jl"))
vis = Visualizer()
open(vis)
# render(vis)

# get model
function get_model(name::String; model_name::String = name, surf::String = "flat", dynamics::String="dynamics")
	#TODO: assert model exists
	path = joinpath(@__DIR__, "..", "dynamics", name)
	# include(joinpath(path, "model.jl"))
	model = eval(Symbol(model_name * (surf != "flat" ? "_" * surf : "")))
	instantiate_base!(model, joinpath(path, dynamics, "base.jld2"))
	instantiate_dynamics!(model, joinpath(path, dynamics, "dynamics.jld2"))
	instantiate_residual!(model, joinpath(path, surf, "residual.jld2"))
	# instantiate_linearized!(model, joinpath(path, surf, "linearized.jld2"))
	@load joinpath(path, surf, "sparse_jacobians.jld2") rz_sp rθ_sp
	model.spa.rz_sp = rz_sp
	model.spa.rθ_sp = rθ_sp
	return model
end
model = get_model("quadrupedlinear")

# get trajectory
# ref_traj = get_trajectory("quadrupedlinear", "quadruped_v2_mirror_gait_fast2", load_type = :split_traj_alt)
# ref_traj = get_trajectory("quadrupedlinear", "quadruped_v2_mirror_gait_fast", load_type = :split_traj_alt)
ref_traj = get_trajectory("quadrupedlinear", "gait1_mit_2.5percent", load_type = :split_traj_alt)
# ref_traj = get_trajectory("quadrupedlinear", "gait0_mit_10percent", load_type = :split_traj_alt)
ref_traj_copy = deepcopy(ref_traj)


function check_traj(model::ContactDynamicsModel, traj::ContactTraj)
    nq = model.dim.q
    nz = num_var(model)
    d = [zeros(nq) for t=1:traj.H]
    r = [zeros(nz) for t=1:traj.H]
    for t = 1:ref_traj.H
        k = kinematics(model, traj.q[t+2])
        λ = contact_forces(model, traj.γ[t], traj.b[t], traj.q[t+2], k)
        d[t] = dynamics(model, traj.h, traj.q[t], traj.q[t+1], traj.u[t], traj.w[t], λ, traj.q[t+2])
        r[t] = residual(model, traj.z[t], traj.θ[t], [1e-8])
    end
    return d, r
end
ref_traj
nz = num_var(model)
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
N_sample = 1
H_mpc = 20
h_sim = h / N_sample
H_sim = 300 #4000 #3000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(H_mpc, model.dim,
    q = [Diagonal(1e-1 * [1e-1; 1e0; 1e0; 1e0*ones(3); 1.0 * ones(model.dim.q-6)]) for t = 1:H_mpc],
	v = [Diagonal(1e-2 * [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0 * ones(model.dim.q-6)]) for t = 1:H_mpc],
    u = [Diagonal(1e-2 * [1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0, 1e0]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])



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

p = linearized_mpc_policy(ref_traj, model, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        solver = :ldl_solver,
        max_iter = 5),
    mpc_opts = LinearizedMPCOptions(
		# live_plotting=true
		))


q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6,
        diff_sol = true),
    sim_opts = ContactControl.SimulatorOptions(
		warmstart = true))

time = @elapsed status = ContactControl.simulate!(sim)
# @elapsed status = ContactControl.simulate!(sim)
# @profiler status = ContactControl.simulate!(sim)

# plot_lines!(vis, model, sim.traj.q[1:25:end])
plot_surface!(vis, model.env, ylims=[0.3, -0.3])
anim = visualize_robot!(vis, model, ref_traj, name=:Ref, sample=1, α=0.5)
anim = visualize_robot!(vis, model, sim.traj, anim=anim, sample=1)
# anim = visualize_robot!(vis, model, sim.traj, anim=anim)
# anim = visualize_force!(vis, model, sim.traj, anim=anim, h=h_sim)

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


filename = "tablebot_working"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)
