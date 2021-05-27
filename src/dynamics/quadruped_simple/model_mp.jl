struct QuadrupedLinear11V2{I, T} <: Model{I, T}
    n::Int
    m::Int
    d::Int

    g
    μ

	mb # body mass
	mf # foot mass
	Ix # body inertia (x)
	Iy # body inertia (y)
	Iz # body inertia (z)

	l_torso
	w_torso

    # joint limits
    qL
    qU

    # torque limits
    uL
    uU

    nq
    nu
    nc
    nf
    nb
    ns

    idx_u
    idx_λ
    idx_b
    idx_ψ
    idx_η
    idx_s

	joint_friction
end

# Dimensions
nq = 3 + 3 + 3 * 4        # configuration dimension
nu = 4 * 3                # control dimension
nc = 4                    # number of contact points
nf = 4                    # number of parameters for friction cone
nb = nc * nf
ns = 1

# World parameters
μ = 0.5      # coefficient of friction
g = 9.81     # gravity

# ~Unitree A1
mb = 1.0
mf = 0.01
Ix = 0.1
Iy = 0.1
Iz = 0.1

l_torso = 0.5
w_torso = 0.25

n = 2 * nq
m = nu + nc + nb + nc + nb + ns
d = 0

idx_u = (1:nu)
idx_λ = nu .+ (1:nc)
idx_b = nu + nc .+ (1:nb)
idx_ψ = nu + nc + nb .+ (1:nc)
idx_η = nu + nc + nb + nc .+ (1:nb)
idx_s = nu + nc + nb + nc + nb .+ (1:ns)

qL = -Inf * ones(nq)
qU = Inf * ones(nq)

uL = -100.0 * ones(nu)
uU = 100.0 * ones(nu)

function lagrangian_derivatives(model, q, v)
	D1L = -1.0 * C_func(model, q, v)
    D2L = M_func(model, q) * v
	return D1L, D2L
end


# Methods
function M_func(model::QuadrupedLinear11V2, q)
	M = Diagonal([(model.mb + 4 * model.mf) * ones(3) ;
				  model.Ix; model.Iy; model.Iz;
				  model.mf * ones(12)])
	return M
end

function C_func(model::QuadrupedLinear11V2, q, q̇)
	SVector{18}([0.0, 0.0, model.g * (model.mb + 4 * model.mf),
				 cross(q̇[4:6], Diagonal([model.Ix, model.Iy, model.Iz]) * q̇[4:6])...,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf,
				 0.0, 0.0, model.g * model.mf])
end

function kinematics(model::QuadrupedLinear11V2, q)
	p_torso = q[1:3]

	rot = MRP(q[4:6]...)

	# # feet positions in body frame
	# pb1 = q[6 .+ (1:3)]
	# pb2 = q[9 .+ (1:3)]
	# pb3 = q[12 .+ (1:3)]
	# pb4 = q[15 .+ (1:3)]
	#
	# # feet positions in world frame
	# pw1 = p_torso + rot * pb1
	# pw2 = p_torso + rot * pb2
	# pw3 = p_torso + rot * pb3
	# pw4 = p_torso + rot * pb4

	pw1 = q[6 .+ (1:3)]
	pw2 = q[9 .+ (1:3)]
	pw3 = q[12 .+ (1:3)]
	pw4 = q[15 .+ (1:3)]

	SVector{12}([pw1; pw2; pw3; pw4])
end

function ϕ_func(model::QuadrupedLinear11V2, q)
	k = kinematics(model, q)

	@SVector [k[3], k[6], k[9], k[12]]
end

function skew(x)
	SMatrix{3,3}([0.0 -x[3] x[2];
	               x[3] 0.0 -x[1];
				   -x[2] x[1] 0.0])
end

function B_func(model::QuadrupedLinear11V2, q)
	p_torso = q[1:3]

	rot = MRP(q[4:6]...)

	# r in world frame
	r1 = q[6 .+ (1:3)] - p_torso
	r2 = q[9 .+ (1:3)] - p_torso
	r3 = q[12 .+ (1:3)] - p_torso
	r4 = q[15 .+ (1:3)] - p_torso

	z3 = zeros(3, 3)

	SMatrix{18, 12}([I I I I;
	                transpose(rot) * skew(r1) transpose(rot) * skew(r2) transpose(rot) * skew(r3) transpose(rot) * skew(r4);
					-I z3 z3 z3;
					z3 -I z3 z3;
					z3 z3 -I z3;
					z3 z3 z3 -I])
end

function N_func(model::QuadrupedLinear11V2, q)
	p_torso = q[1:3]

	rot = MRP(q[4:6]...)

	# r in world frame
	r1 = q[6 .+ (1:3)] - p_torso
	r2 = q[9 .+ (1:3)] - p_torso
	r3 = q[12 .+ (1:3)] - p_torso
	r4 = q[15 .+ (1:3)] - p_torso

	z3 = zeros(3, 3)

	N = [zeros(6, 12);
	     I z3 z3 z3;
	     z3 I z3 z3;
	     z3 z3 I z3;
	     z3 z3 z3 I]

	idx = collect([3, 6, 9, 12])
	transpose(N[:, idx])
end

function _P_func(model::QuadrupedLinear11V2, q)
	p_torso = q[1:3]

	rot = MRP(q[4:6]...)

	# r in world frame
	r1 = q[6 .+ (1:3)] - p_torso
	r2 = q[9 .+ (1:3)] - p_torso
	r3 = q[12 .+ (1:3)] - p_torso
	r4 = q[15 .+ (1:3)] - p_torso

	z3 = zeros(3, 3)

	P = [zeros(6, 12);
		 I z3 z3 z3;
		 z3 I z3 z3;
		 z3 z3 I z3;
		 z3 z3 z3 I]

	idx = collect([1, 2, 4, 5, 7, 8, 10, 11])
	transpose(P[:, idx])
end

function P_func(model::QuadrupedLinear11V2, q)
	_P = _P_func(model, q)

	map = [1.0 0.0;
		   0.0 1.0;
		   -1.0 0.0;
		   0.0 -1.0]

	return [map * _P[1:2, :];
	        map * _P[3:4, :];
			map * _P[5:6, :];
			map * _P[7:8, :]]
end

function friction_cone(model::QuadrupedLinear11V2, u)
	λ = view(u, model.idx_λ)
	b = view(u, model.idx_b)

	return @SVector [model.μ * λ[1] - sum(view(b, 1:4)),
					 model.μ * λ[2] - sum(view(b, 5:8)),
					 model.μ * λ[3] - sum(view(b, 9:12)),
					 model.μ * λ[4] - sum(view(b, 13:16))]
end

function maximum_dissipation(model::QuadrupedLinear11V2{Discrete, FixedTime}, x⁺, u, h)
	q3 = view(x⁺, model.nq .+ (1:model.nq))
	q2 = view(x⁺, 1:model.nq)
	ψ = view(u, model.idx_ψ)
	ψ_stack = [ψ[1] * ones(4); ψ[2] * ones(4); ψ[3] * ones(4); ψ[4] * ones(4)]
	η = view(u, model.idx_η)
	return P_func(model, q3) * (q3 - q2) / h + ψ_stack - η
end

function fd(model::QuadrupedLinear11V2{Discrete, FixedTime}, x⁺, x, u, w, h, t)
	q3 = view(x⁺, model.nq .+ (1:model.nq))
	q2⁺ = view(x⁺, 1:model.nq)
	q2⁻ = view(x, model.nq .+ (1:model.nq))
	q1 = view(x, 1:model.nq)
	u_ctrl = view(u, model.idx_u)
	λ = view(u, model.idx_λ)
	b = view(u, model.idx_b)

	qm1 = 0.5 * (q1 + q2⁺)
    vm1 = (q2⁺ - q1) / h
    qm2 = 0.5 * (q2⁺ + q3)
    vm2 = (q3 - q2⁺) / h


	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	[q2⁺ - q2⁻;
     (0.5 * h * D1L1 + D2L1 + 0.5 * h * D1L2 - D2L2
     + B_func(model, qm2) * SVector{12}(u_ctrl)
     + transpose(N_func(model, q3)) * SVector{4}(λ)
     + transpose(P_func(model, q3)) * SVector{16}(b))]
end

function fd(model::QuadrupedLinear11V2{Discrete, FreeTime}, x⁺, x, u, w, h, t)
	q3 = view(x⁺, model.nq .+ (1:model.nq))
	q2⁺ = view(x⁺, 1:model.nq)
	q2⁻ = view(x, model.nq .+ (1:model.nq))
	q1 = view(x, 1:model.nq)
	u_ctrl = view(u, model.idx_u)
	λ = view(u, model.idx_λ)
	b = view(u, model.idx_b)
	h = u[end]

	qm1 = 0.5 * (q1 + q2⁺)
    vm1 = (q2⁺ - q1) / h
    qm2 = 0.5 * (q2⁺ + q3)
    vm2 = (q3 - q2⁺) / h

	D1L1, D2L1 = lagrangian_derivatives(model, qm1, vm1)
	D1L2, D2L2 = lagrangian_derivatives(model, qm2, vm2)

	[q2⁺ - q2⁻;
     (0.5 * h * D1L1 + D2L1 + 0.5 * h * D1L2 - D2L2
     + B_func(model, qm2) * SVector{12}(u_ctrl)
     + transpose(N_func(model, q3)) * SVector{4}(λ)
     + transpose(P_func(model, q3)) * SVector{16}(b))]
end

function maximum_dissipation(model::QuadrupedLinear11V2{Discrete, FreeTime}, x⁺, u, h)
	q3 = view(x⁺, model.nq .+ (1:model.nq))
	q2 = view(x⁺, 1:model.nq)
	ψ = view(u, model.idx_ψ)
	ψ_stack = [ψ[1] * ones(4); ψ[2] * ones(4); ψ[3] * ones(4); ψ[4] * ones(4)]
	η = view(u, model.idx_η)
	h = u[end]
	return P_func(model, q3) * (q3 - q2) / h + ψ_stack - η
end

model = QuadrupedLinear11V2{Discrete, FixedTime}(n, m, d,
			  g,
			  μ,
			  mb,
			  mf,
			  Ix,
			  Iy,
			  Iz,
			  l_torso,
			  w_torso,
			  qL, qU,
			  uL, uU,
			  nq,
			  nu,
			  nc,
			  nf,
			  nb,
			  ns,
			  idx_u,
			  idx_λ,
			  idx_b,
			  idx_ψ,
			  idx_η,
			  idx_s,
			  nothing)

function visualize!(vis, model::QuadrupedLinear11V2, q;
      r = 0.05, Δt = 0.1)

	default_background!(vis)

	setobject!(vis["torso"],
    	Rect(Vec(-model.l_torso, -model.w_torso, -0.05),Vec(2.0 * model.l_torso, 2.0 * model.w_torso, 0.05)),
    	MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0)))

	feet1 = setobject!(vis["feet1"], Sphere(Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet2 = setobject!(vis["feet2"], Sphere(Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet3 = setobject!(vis["feet3"], Sphere(Point3f0(0),
		convert(Float32, r)),
		MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	feet4 = setobject!(vis["feet4"], Sphere(Point3f0(0),
        convert(Float32, r)),
        MeshPhongMaterial(color = RGBA(1.0, 165.0 / 255.0, 0, 1.0)))

	anim = MeshCat.Animation(convert(Int, floor(1.0 / Δt)))

	T = length(q)
	p_shift = [0.0, 0.0, r]

	for t = 1:T
		MeshCat.atframe(anim, t) do
			rot = MRP(q[t][4:6]...)

			p_torso = [q[t][1]; q[t][2]; q[t][3]] + p_shift
			p_foot1 = q[t][6 .+ (1:3)] + p_shift
			p_foot2 = q[t][9 .+ (1:3)] + p_shift
			p_foot3 = q[t][12 .+ (1:3)] + p_shift
			p_foot4 = q[t][15 .+ (1:3)] + p_shift

			settransform!(vis["torso"], Translation(p_torso))
			settransform!(vis["feet1"], Translation(p_foot1))
			settransform!(vis["feet2"], Translation(p_foot2))
			settransform!(vis["feet3"], Translation(p_foot3))
			settransform!(vis["feet4"], Translation(p_foot4))
		end
	end

	settransform!(vis["/Cameras/default"],
	    compose(Translation(0.0, 0.0, -1.0), LinearMap(RotZ(-pi / 2.0))))

	MeshCat.setanimation!(vis, anim)
end

function initial_configuration(model::QuadrupedLinear11V2)
	q1 = zeros(model.nq)

	# position
	q1[3] = 0.5

	# orientation
	mrp = MRP(RotZ(0.0))
	q1[4:6] = [mrp.x; mrp.y; mrp.z]

	# feet positions (in body frame)
	q1[7:9] = [model.l_torso; model.w_torso; 0.0]
	q1[10:12] = [model.l_torso; -model.w_torso; 0.0]
	q1[13:15] = [-model.l_torso; model.w_torso; 0.0]
	q1[16:18] = [-model.l_torso; -model.w_torso; 0.0]

	return q1
end

# q_init = initial_configuration(model)
#
# include(joinpath(module_dir(), "models/visualize.jl"))
# vis = Visualizer()
# render(vis)
# visualize!(vis, model, [q_init], Δt = 0.1)
