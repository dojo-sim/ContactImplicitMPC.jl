
"""
	instantiate_residual!(model,
		path::AbstractString="model/residual.jld2")
Evaluates the residual expressions to generate functions, stores them into the model.
"""
function instantiate_residual!(fct::ResidualMethods, expr::Dict{Symbol,Expr};
	jacobians = :full)

	fct.r!  = eval(expr[:r])

	if jacobians == :full
		fct.rz! = eval(expr[:rz])
		fct.rθ! = eval(expr[:rθ])
	else
		# fct.rz! = rz_approx!
		# fct.rθ! = rθ_approx!
	end

	return nothing
end


include(joinpath(@__DIR__, "..", "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)



# get hopper model
model = get_model("flamingo")

nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# time
H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 10
h_sim = h / N_sample
H_sim = 800

# barrier parameter
κ_mpc = 1.0e-4

cost = CostFunction(H_mpc, model.dim,
    q = [Diagonal(1e-1 * [1.0, 0.01, 0.05, 1.5, 1.5, .15, .15, .0005, .0005]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [10; 1; 10; ones(nu-5); 10; 10]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-3 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.b)) for t = 1:H_mpc])




# q0_sim = SVector{model.dim.q}([0.0, 0.828, -0.00, -0.4, 0.2, -0.5, 0.9, π/2, π/2])
# q1_sim = SVector{model.dim.q}([0.0, 0.828, -0.00, -0.4, 0.2, -0.5, 0.9, π/2, π/2])
q0_sim = SVector{model.dim.q}([0.0, 0.849, -0.00, 0.1, 0.295, -0.3, 0.1, π/2, π/2])
q1_sim = SVector{model.dim.q}([0.0, 0.849, -0.00, 0.1, 0.295, -0.3, 0.1, π/2, π/2])
visualize_robot!(vis, model, [q0_sim])
kinematics_3(model, q0_sim, body=:foot_1, mode=:toe)[2]
kinematics_3(model, q0_sim, body=:foot_2, mode=:toe)[2]


p = flamingo_policy(model, h_sim, qref=q0_sim)

sim = ContactControl.simulator(model, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = ContactControl.InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = ContactControl.SimulatorOptions(warmstart = true))
model.linearized
@time status = ContactControl.simulate!(sim)





plt = plot(layout=(3,1), legend=false)
# plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
#     color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[1:1] for q in sim.traj.q])...)',
	color=:blue, linewidth=1.0)
# plot!(plt[1,1], hcat(Vector.([q[2:2] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
# plot!(plt[1,1], hcat(Vector.([q[3:3] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
# plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][1:nu], N_sample) for i=1:H]...))...)',
    # color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[1:nu] for u in sim.traj.u]*N_sample)...)',
	color=:blue, linewidth=1.0)
plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)',
	color=:blue, linewidth=1.0)

anim = visualize_robot!(vis, model, sim.traj, name=:mpc)
visualize_force!(vis, model, sim.traj, name=:mpc, anim=anim)

filename = "flamingo_pd"
filename = "flamingo_pd"
filename = "flamingo_pd"
filename = "flamingo_pd"
filename = "flamingo_pd"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)

# const ContactControl = Main

all(p.contact)



function kinematic_map(model::Flamingo, q; body=:foot_1)
	#TODO need to compute wrt to the center of pressure not com
	x, z = q[1:2] - kinematics_3(model, q, body=body, mode=:com)
	# θ = q[3] # torso angle
	# s = [x,z,θ]
	s = [x,z]
	return s
end
function jacobian_map(model::Flamingo, q; body=:foot_1)
	k(q) = kinematic_map(model, q, body=body)
	J = ForwardDiff.jacobian(k, q)[:,3:end]
	return J
end

function virtual_actuator_torque2(model::Flamingo, q::AbstractVector,
		f::AbstractVector; body=:foot_1)
	fx, fz = f
	J = jacobian_map(model, q, body=body)
	if body == :foot_1
		iJ = [2,3,6]
	elseif body == :foot_2
		iJ = [4,5,7]
	else
		@error "incorrect body specification"
	end
	# @show J
	τ = J[:,iJ]'*f
	return τ
end


k = kinematic_map(model, q0_sim, body=:foot_2)
J = jacobian_map(model, q0_sim, body=:foot_2)


f = rand(2)
J = jacobian_map(model, q0_sim, body=:foot_1)
virtual_actuator_torque(model, f, body=:foot_1)
fx, fz = f
if body == :foot_1
	iJ = [2,3,6]
elseif body == :foot_2
	iJ = [4,5,7]
else
	@error "incorrect body specification"
end
τ = J[:,iJ]'*f
