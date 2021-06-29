include(joinpath(@__DIR__, "..", "dynamics", "flamingo", "visuals.jl"))
T = Float64
vis = Visualizer()
# render(vis)
open(vis)

s = get_simulation("flamingo", "flat_2D_lc", "flat")
model = s.model
env = s.env
const ContactControl = Main
ref_traj = deepcopy(ContactControl.get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

H = ref_traj.H
h = ref_traj.h
N_sample = 5
H_mpc = 15
h_sim = h / N_sample
H_sim = 600 # 35000

# barrier parameter
κ_mpc = 1.0e-4

obj = TrackingVelocityObjective(model, env, H_mpc,
    v = [Diagonal(1e-3 * [1e0,1,1e4,1,1,1,1,1e4,1e4]) for t = 1:H_mpc],
    q = [Diagonal(1e-1 * [3e2, 1e-6, 3e2, 1, 1, 1, 1, 0.1, 0.1]) for t = 1:H_mpc],
    u = [Diagonal(3e-1 * [0.1; 0.1; 0.3; 0.3; ones(model.dim.u-6); 2; 2]) for t = 1:H_mpc],
    γ = [Diagonal(1.0e-100 * ones(model.dim.c)) for t = 1:H_mpc],
    b = [Diagonal(1.0e-100 * ones(model.dim.c * friction_dim(env))) for t = 1:H_mpc])

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	ip_type = :interior_point,
	# mode = :configurationforce,
	mode = :configuration,
    n_opts = NewtonOptions(
        r_tol = 3e-4,
        max_iter = 5,
		solver = :lu_solver,),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        altitude_impact_threshold = 0.02,
        altitude_verbose = true,
        )
    )

p = linearized_mpc_policy(ref_traj, s, obj,
    H_mpc = H_mpc,
    N_sample = N_sample,
    κ_mpc = κ_mpc,
	# mode = :configurationforce,
	mode = :configuration,
	ip_type = :mehrotra,
    n_opts = NewtonOptions(
		r_tol = 3e-4,
		max_iter = 5,
		max_time = ref_traj.h, # HARD REAL TIME
		solver = :lu_solver,),
    mpc_opts = LinearizedMPCOptions(
        # live_plotting=true,
        # altitude_update = true,
        altitude_impact_threshold = 0.02,
        altitude_verbose = true,
        ),
	ip_opts = MehrotraOptions(
		max_iter_inner = 100,
		verbose = true,
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

# sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
#     p = p,
#     ip_opts = MehrotraOptions(
#         r_tol = 1.0e-8,
#         κ_init = 1.0e-8,
#         κ_tol = 2.0e-8),
#     sim_opts = SimulatorOptions(warmstart = true),
# 	ip_type = :mehrotra,
#     )
sim = simulator(s, q0_sim, q1_sim, h_sim, H_sim,
    p = p,
    ip_opts = InteriorPointOptions(
        r_tol = 1.0e-8,
        κ_init = 1.0e-8,
        κ_tol = 2.0e-8),
    sim_opts = SimulatorOptions(warmstart = true),
	ip_type = :interior_point,
    )

telap = @elapsed status = simulate!(sim, verbose = true)
# @profiler status = simulate!(sim, verbose = true)
# telap = 2.75 # H_sim = 600
H_sim * h / (telap - 1.10)


l = 9
lu = 1
plt = plot(layout=(3,1), legend=false)
plot!(plt[1,1], hcat(Vector.(vcat([fill(ref_traj.q[i], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[1,1], hcat(Vector.([q[l:l] for q in sim.traj.q])...)', color=:blue, linewidth=1.0)
plot!(plt[2,1], hcat(Vector.(vcat([fill(ref_traj.u[i][lu:lu], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[3,1], hcat(Vector.(vcat([fill(ref_traj.γ[i][1:nc], N_sample) for i=1:H]...))...)',
    color=:red, linewidth=3.0)
plot!(plt[2,1], hcat(Vector.([u[lu:lu] for u in sim.traj.u]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([γ[1:nc] for γ in sim.traj.γ]*N_sample)...)', color=:blue, linewidth=1.0)
# plot!(plt[3,1], hcat(Vector.([b[1:nb] for b in sim.traj.b]*N_sample)...)', color=:red, linewidth=1.0)

plot_surface!(vis, env, xlims=[-0.5, 1.5], ylims = [-0.5, 0.5])
# anim = visualize_robot!(vis, model, sim.traj, sample=10)
anim = visualize_meshrobot!(vis, model, sim.traj, sample=10)
anim = visualize_force!(vis, model, env, sim.traj, anim=anim, h=h_sim, sample=10)


# Test robustness
s = get_simulation("flamingo", "flat_2D_lc", "flat")
ref_traj = deepcopy(get_trajectory(s.model, s.env,
    joinpath(module_dir(), "src/dynamics/flamingo/gaits/gait_forward_36_4.jld2"),
    load_type = :split_traj_alt))

im_traj = ImplicitTraj(ref_traj, s,
	ip_type = :mehrotra,
	opts = MehrotraOptions(
		max_iter_inner = 100,
		κ_init = 1e-8,
		κ_tol = 2.0 * 1e-8,
		r_tol = 1.0e-8,
		diff_sol = true,
		verbose = true,
		solver = :empty_solver))

im_traj.ip[1].z .= ref_traj.z[1]
im_traj.ip[1].θ .= ref_traj.θ[1] .+ 1.0*[ones(length(ref_traj.θ[1])-2); zeros(2)]
interior_point_solve!(im_traj.ip[1])
im_traj.ip[1].iterations

cnt = 0
itl = []
sul = 0
for tt = 1:1000
	t = (tt % ref_traj.H) + 1
    Random.seed!(t)
	im_traj.ip[t].z .= ref_traj.z[t]
	im_traj.ip[t].θ .= ref_traj.θ[t] .+ 1e-3*[0.5 .- rand(length(ref_traj.θ[t])-2); zeros(2)]
	interior_point_solve!(im_traj.ip[t])
	it = im_traj.ip[t].iterations
	su = it < 100
    sul += Int(su)
    push!(itl, it)
end
sul
mean(itl)



filename = "flamingo_mehrotra"
MeshCat.convert_frames_to_video(
    "/home/simon/Downloads/$filename.tar",
    "/home/simon/Documents/$filename.mp4", overwrite=true)

convert_video_to_gif(
    "/home/simon/Documents/$filename.mp4",
    "/home/simon/Documents/$filename.gif", overwrite=true)

































	μ = z_y1'*z_y2 / length(z_y1)
	σ = (μaff / μ)^3




	ip0 = deepcopy(p.im_traj.ip[1])
	centering(ip0.z, ip0.Δaff, ip0.iy1, ip0.iy2, 1.0)
	@code_warntype centering(ip0.z, ip0.Δaff, ip0.iy1, ip0.iy2, 1.0)
	@benchmark centering($ip0.z, $ip0.Δaff, $ip0.iy1, $ip0.iy2, $1.0)


	step_length_11(ip0)
	@code_warntype step_length_11(ip0)
	@benchmark step_length_11($ip0)



	function least_squares!(ip::Mehrotra{T,nx,ny,R,RZ,Rθ}) where {T,nx,ny,R,RZ,Rθ}
		least_squares!(ip.z, ip.θ, ip.r, ip.rz)
		return nothing
	end

	function least_squares!(z::Vector{T}, θ::AbstractVector{T}, r::RLin{T}, rz::RZLin{T}) where {T}
		δθ = r.θ0 - θ
		δrdyn = r.rdyn0 - r.rθdyn * δθ
		δrrst = r.rrst0 - r.rθrst * δθ

		δw1 = rz.A1 * δrdyn + rz.A2 * δrrst
		δw2 = rz.A3 * δrdyn + rz.A4 * δrrst
		δw3 = rz.A5 * δrdyn + rz.A6 * δrrst

		@. @inbounds z[r.ix]  .= r.x0  .+ δw1
		@. @inbounds z[r.iy1] .= r.y10 .+ δw2
		@. @inbounds z[r.iy2] .= r.y20 .+ δw3
		return nothing
	end



	ip0 = deepcopy(p.im_traj.ip[1])
	least_squares!(ip0)
	@code_warntype least_squares!(ip0)
	@benchmark least_squares!($ip0)



	a = 10

	nx = length(r.ix)
	ny = length(r.iy1)

	δrdyn = r.rdyn0 - r.rθdyn * δθ
	δrrst = r.rrst0 - r.rθrst * δθ

	δw1 = rz.A1 * δrdyn + rz.A2 * δrrst
	δw2 = rz.A3 * δrdyn + rz.A4 * δrrst
	δw3 = rz.A5 * δrdyn + rz.A6 * δrrst

	z[ix]  .= r.x0  + δw1
	z[iy1] .= r.y10 + δw2
	z[iy2] .= r.y20 + δw3
	comp && println("**** z+wt:", scn(norm(z), digits=4))

	z .= initial_state!(z, ix, iy1, iy2, comp = comp)

	function test_ip1(ip::Mehrotra{T}, r::RLin) where {T}
		test_1(ip.z, ip.iy1, r.y10)
		return nothing
	end

	function test_1(z::Vector{T}, iy1::SVector{ny,Int}, y10::SVector{ny,T}) where {ny,T}
		z[iy1] = y10
		return nothing
	end


	ip0 = deepcopy(p.im_traj.ip[1])
	test_ip1(ip0, ip0.r)
	@benchmark test_ip1($ip0, $ip0.r)
	@benchmark test_ip1($ip0)
	@code_warntype test_ip1(ip0, ip0.r)




	# nz = num_var(s.model, s.env)
	# nθ = num_data(s.model)
	# clearconsole()
	# Random.seed!(100)
	# ip_ = deepcopy(sim.p.im_traj.ip[2])
	# ip_.opts.max_iter_inner = 10
	#
	# ip_
	#
	# z0_ = zeros(nz)
	# θ0_ = deepcopy(ip_.r.θ0) + 1e-2*[rand(nθ-2); zeros(2)]
	# z0_[[ip_.r.ix; ip_.r.iy1; ip_.r.iy2]] = [ip_.r.x0; ip_.r.y10; ip_.r.y20]
	#
	# interior_point_solve!(ip_, z0_, θ0_)
