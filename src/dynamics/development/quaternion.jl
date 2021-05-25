# Visualization
vis = Visualizer()
open(vis)

# Midpoint
q0 = rand(4)
q0 ./= norm(q0)
q0_ = Rotations.UnitQuaternion(q0)

ϕ = 0.05*rand(3)
qϕ = MRP(ϕ...)

qϕ1 = 1 / sqrt(1 + norm(ϕ)^2) .* [1; ϕ]
qϕ1_ = UnitQuaternion(qϕ1)

q1 = rand(4)
q1 ./= norm(q1)
q1_ = Rotations.UnitQuaternion(q1)

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


build_frame!(vis, name=:frame0, col=1.0)
build_frame!(vis, name=:frame1, col=0.5)

set_frame!(vis, q0_, name=:frame0)
set_frame!(vis, q1_, name=:frame1)
