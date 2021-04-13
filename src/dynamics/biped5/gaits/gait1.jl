using Plots

# Model
include_model("biped")

model = free_time_model(model)

# Visualize
# - Pkg.add any external deps from visualize.jl
include(joinpath(pwd(), "models/visualize.jl"))
vis = Visualizer()
render(vis)

function get_q⁺(x)
	view(x, model.nq .+ (1:model.nq))
end

function ellipse_traj(x_start, x_goal, z, T)
	dist = x_goal - x_start
	a = 0.5 * dist
	b = z

	z̄ = 0.0

	x = range(x_start, stop = x_goal, length = T)

	z = sqrt.(max.(0.0, (b^2) * (1.0 .- ((x .- (x_start + a)).^2.0) / (a^2.0))))

	return x, z
end

# Horizon
T = 61
Tm = 31 #convert(Int, floor(0.5 * T))

# Time step
tf = 1.0
h = tf / (T - 1)

# Configurations
# 1: x pos
# 2: z pos
# 3: torso angle (rel. to downward vertical)
# 4: thigh 1 angle (rel. to downward vertical)
# 5: calf 1 (rel. to downward vertical)
# 6: thigh 2 (rel. to downward vertical)
# 7: calf 2 (rel. to downward vertical)
# 8: foot 1 (rel. to downward vertical)
# 9: foot 2 (rel. to downward vertical)
function initial_configuration(model, θ_torso, θ_thigh_1, θ_leg_1, θ_thigh_2)
    q1 = zeros(model.nq)
    q1[3] = θ_torso
    q1[4] = θ_thigh_1
    q1[5] = θ_leg_1
    z1 = model.l_thigh1 * cos(q1[4]) + model.l_calf1 * cos(q1[5])
    q1[6] = θ_thigh_2
    q1[7] = -1.0 * acos((z1 - model.l_thigh2 * cos(q1[6])) / model.l_calf2)
	# q1[6] = θ_thigh_1
	# q1[7] = θ_leg_1
	q1[2] = z1

    # p1 = kinematics_2(model, q1, body = :calf_1, mode = :ee)
    # p2 = kinematics_2(model, q1, body = :calf_2, mode = :ee)
    # @show stride = abs(p1[1] - p2[1])
	#
    # q1[1] = -0.75#-1.0 * p1[1]
	#
    # qM = copy(q1)
    # qM[4] = q1[6]
    # qM[5] = q1[7]
    # qM[6] = q1[4]
    # qM[7] = q1[5]
    # qM[1] = abs(p2[1])
	#
    # pM_1 = kinematics_2(model, qM, body = :calf_1, mode = :ee)
    # pM_2 = kinematics_2(model, qM, body = :calf_2, mode = :ee)
	#
    # qT = copy(q1)
    # qT[1] = 0.75#q1[1] + 1.0#2 * stride
	#
    # pT_1 = kinematics_2(model, qT, body = :calf_1, mode = :ee)
    # pT_2 = kinematics_2(model, qT, body = :calf_2, mode = :ee)

    return q1
end

# q1 = initial_configuration(model, 0.0, pi / 75.0, 0.0, -pi / 100.0)
q1 = initial_configuration(model, -pi / 100.0, pi / 20.0, -pi / 20.0, -pi / 40.0)
pf1 = kinematics_2(model, q1, body = :calf_1, mode = :ee)
pf2 = kinematics_2(model, q1, body = :calf_2, mode = :ee)

strd = 2 * (pf1 - pf2)[1]

qT = copy(q1)
qT[1] += strd
q_ref = linear_interpolation(q1, qT, T+1)
x_ref = configuration_to_state(q_ref)
visualize!(vis, model, q_ref, Δt = h)

# pt1 = kinematics_3(model, q1, body = :foot_1, mode = :toe)
# ph1 = kinematics_3(model, q1, body = :foot_1, mode = :heel)
#
# pt2 = kinematics_3(model, q1, body = :foot_2, mode = :toe)
# ph2 = kinematics_3(model, q1, body = :foot_2, mode = :heel)

T_fix = 10
zh = 0.05
xf2_el, zf2_el = ellipse_traj(pf2[1], pf2[1] + strd, zh, Tm - T_fix)
xf1_el, zf1_el = ellipse_traj(pf1[1], pf1[1] + strd, zh, Tm - T_fix)

zf2 = [[zf2_el[1] for t = 1:T_fix]..., zf2_el..., [zf2_el[end] for t = 1:Tm-1]...]
xf2 = [[xf2_el[1] for t = 1:T_fix]..., xf2_el..., [xf2_el[end] for t = 1:Tm-1]...]
zf1 = [[zf1_el[1] for t = 1:Tm-1 + T_fix]..., zf1_el...]
xf1 = [[xf1_el[1] for t = 1:Tm-1 + T_fix]..., xf1_el...]

p1_ref = [[xf1[t]; zf1[t]] for t = 1:T]
p2_ref = [[xf2[t]; zf2[t]] for t = 1:T]

plot(hcat(p1_ref...)', legend = :topleft)
plot!(hcat(p2_ref...)')


# Control
# u = (τ1..7, λ1..4, β1..8, ψ1..4, η1...8, s1)
# τ1: torso angle
# τ2: thigh 1 angle
# τ3: calf 1
# τ4: thigh 2
# τ5: calf 2
# τ6: foot 1
# τ7: foot 2

# ul <= u <= uu
# u1 = initial_torque(model, q1, h)[model.idx_u] # gravity compensation for current q
_uu = Inf * ones(model.m)
_uu[model.idx_u] .= Inf
_uu[end] = 2.0 * h
_ul = zeros(model.m)
_ul[model.idx_u] .= -Inf
_ul[end] = 0.25 * h
ul, uu = control_bounds(model, T, _ul, _uu)

qL = [-Inf; -Inf; q1[3] - pi / 50.0; q1[4:end] .- pi / 6.0; -Inf; -Inf; q1[3] - pi / 50.0; q1[4:end] .- pi / 6.0]
qU = [Inf; q1[2] + 0.001; 0.0; q1[4:end] .+ pi / 6.0; Inf; Inf; q1[3:end] .+ pi / 6.0]

xl, xu = state_bounds(model, T,
    qL, qU,
    x1 = [Inf * ones(model.nq); q1],
	xT = [Inf * ones(model.nq); qT[1]; Inf * ones(model.nq - 1)])

# Objective
include_objective(["velocity", "nonlinear_stage", "control_velocity"])

x0 = configuration_to_state(q_ref)

# penalty on slack variable
obj_penalty = PenaltyObjective(1.0e4, model.m-1)

# quadratic tracking objective
# Σ (x - xref)' Q (x - x_ref) + (u - u_ref)' R (u - u_ref)
q_penalty = 0.0 * ones(model.nq)
# q_penalty[1] = 1.0
# q_penalty[1] = 1.0
# q_penalty[2] = 1.0
# q_penalty[3] = 10.0
# q_penalty[8] = 10.0
# q_penalty[9] = 10.0
x_penalty = [q_penalty; q_penalty]
obj_control = quadratic_time_tracking_objective(
    [Diagonal(x_penalty) for t = 1:T],
    [Diagonal([1.0e-1 * ones(model.nu)..., 1.0e-1 * ones(model.nc)..., 1.0e-1 * ones(model.nb)..., 1.0e-8 * ones(model.m - model.nu - model.nc - model.nb - 1)..., 0.0]) for t = 1:T-1],
    [x_ref[end] for t = 1:T],
    [[zeros(model.nu); zeros(model.m - model.nu)] for t = 1:T-1],
	1.0)

obj_ctrl_vel = control_velocity_objective(Diagonal([1.0e-1 * ones(model.nu); 1.0e-2 * ones(model.nc + model.nb); zeros(model.m - model.nu - model.nc - model.nb)]))

# quadratic velocity penalty
q_v = 1.0e-2 * ones(model.nq)
# q_v[3] = 100.0
# q_v[3:7] .= 1.0e-3
# q_v[8:9] .= 100.0
obj_velocity = velocity_objective(
    [Diagonal(q_v) for t = 1:T-1],
    model.nq,
    h = h)

function l_foot_height(x, u, t)
	J = 0.0
	q1 = view(x, 1:7)
	q2 = view(x, 7 .+ (1:7))

	if t >= Tm
		pq1 = kinematics_2(model, q1, body = :calf_1, mode = :ee)
		pq2 = kinematics_2(model, q1, body = :calf_1, mode = :ee)
		v = (pq2 - pq1) ./ h
		p_avg = 0.5 * (pq1[1] + pq2[1])
		p_max = max(pq1[1], pq2[1])
		J += 1000.0 * sum((p1_ref[t] - kinematics_2(model, q1, body = :calf_1, mode = :ee)).^2.0)
		J += 1000.0 * sum((p1_ref[t] - kinematics_2(model, q1, body = :calf_1, mode = :ee)).^2.0)
		# J += 1000.0 * v' * v
		# J += 1000.0 * (x_ref[end][1] - q1[1])^2.0
	end

	if t <= Tm
		pq1 = kinematics_2(model, q1, body = :calf_2, mode = :ee)
		pq2 = kinematics_2(model, q2, body = :calf_2, mode = :ee)
		v = (pq2 - pq1) ./ h
		p_avg = 0.5 * (pq1[1] + pq2[1])
		p_max = max(pq1[1], pq2[1])

		J += 1000.0 * sum((p2_ref[t] - kinematics_2(model, q1, body = :calf_2, mode = :ee)).^2.0)
		J += 1000.0 * sum((p2_ref[t] - kinematics_2(model, q2, body = :calf_2, mode = :ee)).^2.0)
		# J += 1000.0 * v' * v
		# J += 1000.0 * (x_ref[end][1] - q1[1])^2.0
	end
	#
	# if t < T
	# 	J += 1.0e-3 * ((q2[1] - q1[1]) / max(1.0e-6, u[model.m]) - 1.0)^2.0
	# end

	return J
end

l_foot_height(x) = l_foot_height(x, nothing, T)

obj_foot_height = nonlinear_stage_objective(l_foot_height, l_foot_height)

obj = MultiObjective([obj_penalty,
                      obj_control,
                      obj_velocity,
					  obj_ctrl_vel,
					  obj_foot_height])
# Constraints
include_constraints(["contact", "loop", "free_time", "stage"])

function pinned1!(c, x, u, t)
    q1 = view(x, 1:7)
	# q2 = view(x, 9 .+ (1:9))
    c[1:2] = p1_ref[t] - kinematics_2(model, q1, body = :calf_1, mode = :ee)
	# c[3:4] = p2_ref[t] - kinematics_3(model, q1, body = :foot_2, mode = :com)
	nothing
end

function pinned2!(c, x, u, t)
    q1 = view(x, 1:7)
    # c[1:2] = p1_ref[t] - kinematics_3(model, q, body = :foot_1, mode = :com)
	c[1:2] = p2_ref[t] - kinematics_2(model, q1, body = :calf_2, mode = :ee)
	nothing
end

# cc = zeros(4)
# xx = zeros(model.n)
# uuu = zeros(model.m)
# pinned!(cc, xx, uuu, 1)
# p1_ref[1]
n_stage = 2
t_idx1 = vcat([t for t = 1:Tm+T_fix]..., T)
t_idx2 = vcat([1:T_fix + 1]..., [t for t = Tm:T]...)

con_pinned1 = stage_constraints(pinned1!, n_stage, (1:0), t_idx1)
con_pinned2 = stage_constraints(pinned2!, n_stage, (1:0), t_idx2)

# constraints!(zeros(con_pinned.n), zeros(prob.num_var), con_pinned, model, prob.prob.idx, h, T)
con_loop = loop_constraints(model, collect([(2:7)...,(9:14)...]), 1, T)
con_contact = contact_constraints(model, T)
con_free_time = free_time_constraints(T)
con = multiple_constraints([con_contact,
	con_free_time,
	con_loop,
	con_pinned1,
	con_pinned2])

# con = multiple_constraints([con_free_time, con_contact, con_pinned1])

# Problem
prob = trajectory_optimization_problem(model,
               obj,
               T,
               xl = xl,
               xu = xu,
               ul = ul,
               uu = uu,
               con = con)

# trajectory initialization
u0 = [[1.0e-2 * randn(model.nu); 0.01 * randn(model.m - model.nu - 1); h] for t = 1:T-1] # random controls


# Pack trajectories into vector
z0 = pack(x0, u0, prob)

# Solve
# NOTE: run multiple times to get good trajectory
# include_snopt()
@time z̄, info = solve(prob, copy(z0),
    nlp = :ipopt,
	max_iter = 2000,
    tol = 1.0e-3, c_tol = 1.0e-3, mapl = 5,
    time_limit = 60 * 3)

# @time z̄, info = solve(prob, copy(z̄ .+ 0.01 * randn(prob.num_var)),
#     nlp = :SNOPT7,
# 	max_iter = 1000,
#     tol = 1.0e-3, c_tol = 1.0e-3, mapl = 5,
#     time_limit = 60 * 3)
# @load joinpath(@__DIR__, "walker_gait_3.jld2") z̄ x̄ ū q̄ τ̄ λ̄ b̄ h̄

@show check_slack(z̄, prob)
x̄, ū = unpack(z̄, prob)
_tf, _t, h̄ = get_time(ū)
q = state_to_configuration(x̄)
u = [u[model.idx_u] for u in ū]
γ = [u[model.idx_λ] for u in ū]
b = [u[model.idx_b] for u in ū]
ψ = [u[model.idx_ψ] for u in ū]
η = [u[model.idx_η] for u in ū]
h̄ = mean(h̄)

t = 1
# ψ[t] .* friction_cone(model, γ[t], b[t])
# model.idx_ψ
# γ[t]
# b[t]
# ψ[t]
# friction_cone(model, γ[t], b[t])

norm(hcat([γ[t] .* ϕ_func(model, q[t+2]) for t = 1:T-1]...),Inf)
norm(hcat([b[t] .* η[t] for t = 1:T-1]...), Inf)
function friction_cone(model::Biped, λ, b)
	# @show model.μ
	return @SVector [model.μ * λ[1] - sum(b[1:2]),
					 model.μ * λ[2] - sum(b[3:4])]
end
norm(hcat([ψ[t] .* friction_cone(model, γ[t], b[t]) for t = 1:T-1]...), Inf)


@save joinpath(@__DIR__, "biped_gait_alt.jld2") z̄ x̄ ū h̄ q u γ b
# @load joinpath(@__DIR__, "walker_gait.jld2") z̄ x̄ ū q̄ τ̄ λ̄ b̄ h̄

maximum([norm(fd(model, x̄[t+1], x̄[t], ū[t], zeros(model.d), h̄, t)) for t = 1:T-1])
(q[end][1] - q[1][1]) / _tf
vis = Visualizer()
render(vis)
visualize!(vis, model,
	[[x̄[1][1:model.nq] for i = 1:10]...,
	 state_to_configuration(x̄)...,
	 [x̄[end][model.nq .+ (1:model.nq)] for i = 1:10]...], Δt = h̄[1])
# visualize!(vis, model,
#  	[q̄[T+1]], Δt = h̄[1])
# #

plot(hcat(u...)', linetype = :steppost)
plot(hcat([γt[1:2] for γt in γ]...)', linetype = :steppost)

plot!(hcat([bt[1:2] for bt in b]...)', linetype = :steppost)
plot!(hcat([bt[3:4] for bt in b]...)', linetype = :steppost)

_pf1 = [kinematics_2(model, qt, body = :calf_1, mode = :ee) for qt in q]
_pf2 = [kinematics_2(model, qt, body = :calf_2, mode = :ee) for qt in q]

# using Plots

plot(hcat(_pf1...)', title = "foot 1",
 	legend = :topleft, label = ["x" "z"], color = :black, width = 2.0)
plot!(hcat(p1_ref...)', legend = :topleft, color = :red, width = 1.0)

plot(hcat(p2_ref...)', title = "foot 2",
	legend = :topleft, label = ["x" "z"], color = :black, width = 2.0)
plot!(hcat(_pf2...)', title = "foot 2", legend = :bottomright,
	color = :red, width = 1.0)

# # # # @save joinpath(pwd(), "examples/trajectories/walker_steps.jld2") z̄
# # # plot(hcat(ū...)[1:model.nu, :]', linetype = :steppost)
function get_q_viz(q̄; N = 4)
	q_viz = [q̄...]
	shift_vec = zeros(model.nq)
	shift_vec[1] = q̄[end][1]
	for i = 1:N
		# println(shift)
		# shift_vec[1] = strd
		#
		q_update = [q + shift_vec for q in q̄[2:end]]
		push!(q_viz, q_update...)
		shift_vec[1] = q_update[end][1]
	end

	return q_viz
end

q_viz = get_q_viz(q)
render(vis)
visualize!(vis, model,
	q_viz,
	Δt = h̄[1])

# function traj_concat(q̄, ū; N = 3)
# 	u_viz = [ū...]
# 	q_viz = [q̄...]
#
# 	shift_vec = zeros(model.nq)
# 	shift_vec[1] = q̄[end][1]
#
# 	for i = 1:N
# 		push!(u_viz, ū...)
# 		# println(shift)
# 		# shift_vec[1] = strd
# 		#
# 		q_update = [q + shift_vec for q in q̄[3:end]]
# 		push!(q_viz, q_update...)
# 		shift_vec[1] = q_update[end][1]
# 	end
#
# 	return q_viz, u_viz
# end
#
# q_cat, u_cat = traj_concat(q̄, ū, N = 25)
# plot(hcat(q_cat...)', label = "")
# plot(hcat(u_cat...)[1:model.nu, :]', label = "")
# x_cat = configuration_to_state(q_cat)
# u_cat
#
# qp = ones(model.nq)
# qp[1] = 1.0
# qp[2] = 1.0
# Q = [Diagonal([qp; qp]) for t = 1:length(x_cat)]
# R = [Diagonal(ones(model.m)) for t = 1:length(u_cat)]
#
# K, P = tvlqr(model, x_cat, u_cat, nothing, Q, R)
#
# plot(hcat([vec(k) for k in K]...)', label = "")
#
# K_walker = K[1:T-1]
# x̄_walker = x̄[1:T]
# ū_walker = ū[1:T-1]
# h̄_walker = h̄[1:T-1]
#
# @save joinpath(@__DIR__, "walker_lqr.jld2") K_walker x̄_walker ū_walker h̄_walker
# # (q̄[end][1] - q̄[1][1]) / tf_
# #
# # h̄[1]
# #
# model_sim = Walker{Discrete, FixedTime}(n, m, d,
# 			  g_world, μ_world,
# 			  l_torso, d_torso, m_torso, J_torso,
# 			  l_thigh, d_thigh, m_thigh, J_thigh,
# 			  l_calf, d_calf, m_calf, J_calf,
# 			  l_foot, d_foot, m_foot, J_foot,
# 			  l_thigh, d_thigh, m_thigh, J_thigh,
# 			  l_calf, d_calf, m_calf, J_calf,
# 			  l_foot, d_foot, m_foot, J_foot,
# 			  qL, qU,
# 			  uL, uU,
# 			  nq,
# 			  nu,
# 			  nc,
# 			  nf,
# 			  nb,
# 			  ns,
# 			  idx_u,
# 			  idx_λ,
# 			  idx_b,
# 			  idx_ψ,
# 			  idx_η,
# 			  idx_s,
# 			  joint_friction)
#
# x_sim = [x̄[1]]
# q_sim = [x̄[1][1:model.nq], x̄[2][model.nq .+ (1:model.nq)]]
# include(joinpath(pwd(), "src/contact_simulator/simulator.jl"))
#
# for t = 1:T-1
# 	_x = step_contact(model_sim, x_sim[end], ū[t][1:model.nu], zeros(model.d), h̄,
# 	        tol_c = 1.0e-5, tol_opt = 1.0e-5, tol_s = 1.0e-4, nlp = :ipopt)
# 	push!(x_sim, _x)
# 	push!(q_sim, x_sim[end][model.nq .+ (1:model.nq)])
# end
#
# plot(hcat(q...)[:, 1:length(q_sim)]', color = :black, width = 2.0)
# plot!(hcat(q_sim...)', color = :red, width = 1.0)
# vis = Visualizer()
# render(vis)
# # open(vis)
# visualize!(vis, model, q_sim, Δt = h̄[1])
