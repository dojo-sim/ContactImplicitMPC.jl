model = get_model("particle")
nq = model.dim.q
nu = model.dim.u
nw = model.dim.w
nc = model.dim.c
nb = model.dim.b

# Visualizer
vis = Visualizer()
open(vis)

function visualize!(vis, model::Particle, q;
	Δt = 0.1, r = 0.25)

	default_background!(vis)

    setobject!(vis["particle"],
		# GeometryBasics.Sphere(GeometryBasics.Point3f0(0),
		# convert(Float32, r)),
		MeshCat.HyperRectangle(MeshCat.Vec(0,0,0),MeshCat.Vec(r,r,r)),
		MeshPhongMaterial(color = RGBA(165.0 / 255.0, 0, 1, 1.0)))

    anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

    for t = 1:length(q)
        MeshCat.atframe(anim, t) do
            settransform!(vis["particle"], MeshCat.Translation(q[t][1:3]...))
        end
    end

	# settransform!(vis["/Cameras/default"],
	#     compose(Translation(-2.5, 7.5, 1.0),LinearMap(RotZ(0.0))))

    MeshCat.setanimation!(vis, anim)
end
visualize!(vis, model, sim.traj.q)



# Demo simulation of particle
H = 20
h = 0.05
u = [+0.2ones(nu) for t=1:H]
w = [-0.05ones(nu) for t=1:H]
q0 = SVector{nq}([0,0,1.])
q1 = SVector{nq}([0,0,1.])
sim = simulator(model, q0, q1, h, H,
    p=open_loop_policy(u, h),
	d=open_loop_disturbances(w, h),
    ip_opts=InteriorPointOptions(r_tol=1e-8, κ_tol=2e-7, κ_init=1e-3),
    sim_opts=SimulatorOptions())

simulate!(sim)
sim
plot(hcat(Vector.(sim.traj.q)...)')

visualize!(vis, model, sim.traj.q)


# EKF
nx = 2nq + nc + nb

Q = zeros(ns,ns) + I(ns)
R = zeros(nq,nq) + I(nq)

P0 = zeros(ns,ns) + I(ns)
q10 = [0,0,1.]
q20 = [0,0,1.]
γ10 = [0.]
b10 = [0, 0, 0, 0.]
x0 = [q10; q20; γ10; b10]
u1 = 0.2*ones(nu)

function meas(x)
	q1, q2, γ1, b1 = unpack_x(x)
	y = q2
	return y
end

function unpack_x(x)
	off = 0
	q1 = x[off .+ (1:nq)]; off += nq
	q2 = x[off .+ (1:nq)]; off += nq
	γ1 = x[off .+ (1:nc)]; off += nc
	b1 = x[off .+ (1:nb)]; off += nb
	return q1, q2, γ1, b1
end

function ekf_step(x, P, u, Q, R)
	# Predict
	# Q = fct(W, df/dw) #TODO
	x̄ = process(x, u,)
	F = df/dx|x̄,u
	H = dh/dx|x̄
	P̄ = F*P*F' + Q #done
	# Update
	y = meas(x̄) + sqrt(R)*randn(nq) #done
	ỹ = y - meas(x̄) #done
	S = H*P̄*H' + R  #done
	K = P̄*H'*inv(S) #done
	x̂ = x̄ + K*ỹ     #done
	P̂ = (I - K*H)*P̄ #done
	return x̂, P̂
end


ekf_step(x0, P0, u1, Q, R)


function process(x, u; w=zeros(nw))
	q0, q1, γ0, b0 = unpack_x(x)
	sim = simulator(model, SVector{nq}(q0), SVector{nq}(q1), h, 1,
		p=open_loop_policy([u], h),
		d=open_loop_disturbances([w], h),
		ip_opts=InteriorPointOptions(r_tol=1e-8, κ_tol=2e-7, κ_init=1e-3, diff_sol=true),
		sim_opts=SimulatorOptions())
	simulate!(sim)
	q2 = sim.traj.q[3]
	γ1 = sim.traj.γ[1]
	b1 = sim.traj.b[1]
	x̄ = [q1; q2; γ1; b1]
	return x̄, sim
end





x = rand(nx)
u = rand(nu)
w = rand(nw)

process(x, u; w=w)
