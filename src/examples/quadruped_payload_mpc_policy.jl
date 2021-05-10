const ContactControl = Main
include(joinpath(@__DIR__, "..", "dynamics", "quadruped", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
render(vis)

# get hopper model
model_payload = get_model("quadruped", surf="payload", dynamics="dynamics_payload")
model_no_payload = get_model("quadruped", surf="flat")

model = get_model("quadruped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd
M_fast(model_payload, zeros(nq))
M_fast(model_no_payload, zeros(nq))

# get trajectory
ref_traj = get_trajectory("quadruped", "gait2", load_type=:split_traj_alt, model=model)
ref_traj_copy = deepcopy(ref_traj)

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 2100 #220 #5000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingObjective(H_mpc, model.dim,
    q = [Diagonal(1e-2 * [10; 0.02; 0.25; 0.25 * ones(nq-3)]) for t = 1:H_mpc],
    u = [Diagonal(3e-2 * ones(model.dim.u)) for t = 1:H_mpc],
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
        # live_plotting=true,
        altitude_update = true,
        altitude_impact_threshold = 0.05,
        altitude_verbose = true,
        )
    )

# Test open loop policy to see the impact of the payload
# p = open_loop_policy(deepcopy(ref_traj.u), N_sample=N_sample)

q1_ref = copy(ref_traj.q[2])
q0_ref = copy(ref_traj.q[1])
q1_sim = SVector{model.dim.q}(q1_ref)
q0_sim = SVector{model.dim.q}(copy(q1_sim - (q1_ref - q0_ref) / N_sample))
@assert norm((q1_sim - q0_sim) / h_sim - (q1_ref - q0_ref) / h) < 1.0e-8

sim_no_payload = simulator(model_no_payload, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = SimulatorOptions(warmstart = true)
    )
@time status = simulate!(sim_no_payload)

sim_payload = simulator(model_payload, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-6,
        κ_tol = 2.0e-6),
    sim_opts = SimulatorOptions(warmstart = true)
    )
@time status = simulate!(sim_payload)

# plt = plot(layout=(3,1), legend=false)
# plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[1,1], hcat(Vector.(sim.traj.q)...)', color=:blue, linewidth=1.0)
# plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
# plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)

plot_surface!(vis, model.env, ylims=[0.3, -0.05])
plot_lines!(vis, model, sim_no_payload.traj.q[1:1:end], name=:NoPayload)
plot_lines!(vis, model, sim_payload.traj.q[1:1:end], name=:Payload)
ext_ref_traj = repeat_ref_traj(ref_traj, model, 7; idx_shift = (1:1))
plot_lines!(vis, model, ext_ref_traj.q, offset=0.025, name=:Ref, col=false)

anim = visualize_meshrobot!(vis, model, sim_no_payload.traj, sample=5, name=:NoPayload)
anim = visualize_meshrobot!(vis, model, sim_payload.traj, anim=anim, sample=5, name=:Payload)
anim = visualize_payload!(vis, model, sim_payload.traj, anim=anim, sample=5, name=:Payload)
anim = visualize_force!(vis, model, sim_no_payload.traj, anim=anim, sample=5, h=h_sim, name=:NoPayload)
anim = visualize_force!(vis, model, sim_payload.traj, anim=anim, sample=5, h=h_sim, name=:Payload)

# Display ghosts
t_ghosts = [1]
mvis_ghosts = []
for (i,t) in enumerate(t_ghosts)
    α = i/(length(t_ghosts)+1)
    name = Symbol("ghost$i")
    mvis = build_meshrobot!(vis, model_payload, name=name, α=α)
    build_payload!(vis, model_payload, name=name, α=α)
    push!(mvis_ghosts, mvis)
end

for (i,t) in enumerate(t_ghosts)
    name = Symbol("ghost$i")
    set_meshrobot!(vis, mvis_ghosts[i], model, sim_payload.traj.q[t], name=name)
    set_payload!(vis, model, sim_payload.traj.q[t], name=name)
end



filename = "quadruped_3kg_vs_ref"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)


function ModifiedMeshFileObject(obj_path::String, material_path::String;
        scale::T = 0.1) where {T}
    obj = MeshFileObject(obj_path)
    rescaled_contents = rescale_contents(obj_path, scale = scale)
    material = select_material(material_path)
    mod_obj = MeshFileObject(
        rescaled_contents,
        obj.format,
        material,
        obj.resources,
        )
    return mod_obj
end

function rescale_contents(obj_path::String; scale::T = 0.1) where T
    lines = readlines(obj_path)
    rescaled_lines = copy(lines)
    for (k,line) in enumerate(lines)
        if length(line) >= 2
            if line[1] == 'v'
                stringvec = split(line, " ")
                vals = map(x -> parse(Float64, x), stringvec[2:end])
                rescaled_vals = vals .* scale
                rescaled_lines[k] = join([stringvec[1]; string.(rescaled_vals)], " ")
            end
        end
    end
    rescaled_contents = join(rescaled_lines, "\r\n")
    return rescaled_contents
end

function select_material(material_path::String)
    mtl_file = open(material_path)
    mtl = read(mtl_file, String)
    return mtl
end

obj_box = joinpath(pwd(), "src/dynamics/quadruped/box/Box.obj")
mtl_box = joinpath(pwd(), "src/dynamics/quadruped/box/Box.mtl")
ctm_box = ModifiedMeshFileObject(obj_box, mtl_box, scale=1.0)

vis = Visualizer()
open(vis)
setobject!(vis["box"],ctm_box)
settransform!(vis["box"], compose(Translation(0.0,0.1,0.1),LinearMap(RotZ(0.0)*RotX(0.0))))

build_payload!(vis, model_payload, name = :Quadruped, r=0.0205, rp=0.13, α=1.0)


tform = LinearMap(0.10*I(3))
settransform!(vis, tform)
