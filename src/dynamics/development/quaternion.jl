# Visualization
vis = Visualizer()
open(vis)

function build_frame!(vis::Visualizer; name::Symbol=:Frame, col::T=1.0, α::T=1.0) where {T}
	matx = MeshPhongMaterial(color = RGBA(col*1.0, 0.0, 0.0, α))
	maty = MeshPhongMaterial(color = RGBA(col*0.9, 0.0, 0.0, α))
	matz = MeshPhongMaterial(color = RGBA(col*0.8, 0.0, 0.0, α))
	setobject!(vis[name][:x], GeometryBasics.HyperRectangle(Vec(0, 0, 0),Vec(1.0, 0.1, 0.1)), matx)
	setobject!(vis[name][:y], GeometryBasics.HyperRectangle(Vec(0, 0, 0),Vec(0.1, 1.0, 0.1)), maty)
    setobject!(vis[name][:z], GeometryBasics.HyperRectangle(Vec(0, 0, 0),Vec(0.1, 0.1, 1.0)), matz)
    return nothing
end

function set_frame!(vis::Visualizer, q::UnitQuaternion; name::Symbol=:Frame)
	settransform!(vis[name], LinearMap(q))
	return
end

# Midpoint
# Reference quaternion
q0 = rand(4)
q0 ./= norm(q0)

# Small rotation
ϕ = 0.15*rand(3)
qϕ = 1 / sqrt(1 + norm(ϕ)^2) .* [1; ϕ]

# Use sqrt form Quaternion Pkg
Qϕ = Quaternion(qϕ...)
Qϕ_mid = sqrt(Qϕ)
qϕ_mid = [Qϕ_mid.s, Qϕ_mid.v1, Qϕ_mid.v2, Qϕ_mid.v3]

qϕ_mid - sqrt_quat(qϕ)
q3 = [1,0,0,0]
q3 - sqrt_quat(q3)



# Full rotation
q1 = L_multiply(qϕ) * q0
# Midpoint rotation
q0_mid = L_multiply(qϕ_mid) * q0
# Check that 1/2 + 1/2 = Full rotations
q1_mid = L_multiply(qϕ_mid) * q0_mid
norm(q1_mid - q1) < 1e-10

build_frame!(vis, name=:frame0, col=1.0)
build_frame!(vis, name=:frame0_mid, col=0.2)
build_frame!(vis, name=:frame1, col=0.5)
build_frame!(vis, name=:frame1_mid, col=0.2)

set_frame!(vis, UnitQuaternion(q0), name=:frame0)
set_frame!(vis, UnitQuaternion(q0_mid), name=:frame0_mid)
set_frame!(vis, UnitQuaternion(q1), name=:frame1)
set_frame!(vis, UnitQuaternion(q1_mid), name=:frame1_mid)



q1L = L_multiply(qϕ) * q0
q1R = R_multiply(q0) * qϕ
norm(q1R - q1L) < 1e-10
(norm(q1L) - 1) < 1e-10
(norm(q1R) - 1) < 1e-10

@variables q[1:4]
sq = sqrt_quat(q)
dsq = Symbolics.jacobian(sq, q)
expr = build_function(dsq, q)[1]
fct = eval(expr)
q_ = rand(4)
@benchmark fct(q_)

@variables q0[1:4]
@variables q1[1:4]
qmid = midpoint(q0, q1)
dqmid = Symbolics.jacobian(qmid, q0)
expr = build_function(dqmid, [q0; q1])[1]
expr = build_function(qmid, [q0; q1])[1]
fct = eval(expr)
q0_ = rand(4)
q1_ = rand(4)
@benchmark fct([q0_; q1_])
