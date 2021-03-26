model = get_model("particle")
nq = model.dim.q
nu = model.dim.u
nw = model.dim.w
nc = model.dim.c
nb = model.dim.b

# Visualizer
vis = Visualizer()
open(vis)


# Demo simulation of particle
Random.seed!(100)
H = 25
h = 0.05
u = [[0.3, 0.9, 0.0] for t=1:H]
w = [-0.15rand(nu) for t=1:H]
q0 = SVector{nq}([0,0,0.5])
q1 = SVector{nq}([0,0,0.5])
sim = simulator(model, q0, q1, h, H,
    p=open_loop_policy(u, h),
	d=open_loop_disturbances(w, h),
    ip_opts=InteriorPointOptions(r_tol=1e-8, κ_tol=2e-7, κ_init=1e-3),
    sim_opts=SimulatorOptions())

simulate!(sim)
plot(hcat(Vector.(sim.traj.q)...)')
visualize!(vis, model, sim.traj.q)


# EKF
nx = 2nq + nc + nb

Q = zeros(nx,nx) + (0.05^2)*I(nx)
R = zeros(nq,nq) + (0.15^2)*I(nq)

P0 = zeros(nx,nx) + (0.01)^2*I(nx)
q0 = SVector{nq}([0,0,0.5])
q1 = SVector{nq}([0,0,0.5])
γ0 = [0.]
b0 = [0, 0, 0, 0.]
x0 = [q0; q1; γ0; b0]
u1 = 0.2*ones(nu)

function ekf_step(x, P, u, y, Q, R)
	# Predict
	# Q = fct(W, df/dw) #TODO
	x̄, F, H = process(x, u) #done
	P̄ = F*P*F' + Q #done
	# Update
	ỹ = y - meas(x̄) #done
	S = H*P̄*H' + R  #done
	K = P̄*H'*inv(S) #done
	x̂ = x̄ + K*ỹ     #done
	P̂ = (I - K*H)*P̄ #done
	return x̂, P̂
end

# x1, P1 = ekf_step(x0, P0, u1, Q, R)
# q1_, q2_, γ1_, b1_ = unpack_x(x1)
# q1_
# q2_
# γ1_
# b1_
#
# plot(Gray.(P1/maximum(P1)))

function ekf_rollout(H::Int, x0, P0, u, w, Q, R)
	Random.seed!(100)
	X = x0 #+ sqrt(P0)*randn(nx) # real state
	X_ = [X]
	x_ = [x0]
	P_ = [P0]
	for t = 1:H
		X, = process(X, u[t], w=w[t])
		y = meas(X) + sqrt(R)*randn(nq) #done
		x̂, P̂ = ekf_step(x_[end], P_[end], u[t], y, Q, R)
		push!(x_, copy(x̂))
		push!(P_, copy(P̂))
		push!(X_, copy(X))
	end
	return X_, x_, P_
end

X_, x_, P_ = ekf_rollout(H, x0, P0, u, w, Q, R)

visualize!(vis, model, sim.traj.q, name=:real_noisy)
visualize!(vis, model, [x[nq .+ (1:nq)] for x in X_], name=:real_noisy2)
visualize_uncertainty!(vis, model,
	[x[nq .+ (1:nq)] for x in x_],
	[P[nq .+ (1:nq), nq .+ (1:nq)] for P in P_], name=:estim)



# Constrained EKF



# filename = "nom_particle"
# MeshCat.convert_frames_to_video(
#     "/home/simon/Downloads/$filename.tar",
#     "/home/simon/Documents/$filename.mp4", overwrite=true)
#
# convert_video_to_gif(
#     "/home/simon/Documents/$filename.mp4",
#     "/home/simon/Documents/$filename.gif", overwrite=true)
