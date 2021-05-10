include(joinpath(@__DIR__, "visuals.jl"))
T = Float64
vis = Visualizer()
open(vis)
# render(vis)

# get hopper model
# model_sim = get_model("quadruped", surf="sinusoidal")
model = get_model("quadruped", surf="flat")
nq = model.dim.q
nu = model.dim.u
nc = model.dim.c
nb = model.dim.b
nd = nq + nc + nb
nr = nq + nu + nc + nb + nd

# get trajectory
ref_traj = get_trajectory("quadruped", "gait1", load_type=:split_traj_alt, model=model)


plot_surface!(vis, model.env, ylims=[-0.5, 0.3])
anim = visualize_meshrobot!(vis, model, ref_traj.q)
visualize_robot!(vis, model, ref_traj.q, anim=anim)

# Test visualizer
plot_surface!(vis, model.env, ylims=[-0.5, 0.3])
build_robot!(vis, model, name=:quad0f)
build_robot!(vis, model, name=:quad0b)
mvis = build_meshrobot!(vis, model, name=:shadow_1, α=0.5)

t = 20
q = ref_traj.q[t] + [0,0,pi/1,0,0,0,0,0,0,0,0]
set_robot!(vis, model, q, name=:quad0f, offset=0.00)
set_robot!(vis, model, q, name=:quad0b, offset=0.264)
set_meshrobot!(vis, mvis, model, q, name=:shadow_1)



vis2 = Visualizer()
open(vis2)

matx = MeshPhongMaterial(color = RGBA(0.0, 0.0, 0.0, 1.0))
maty = MeshPhongMaterial(color = RGBA(0.5, 0.5, 0.5, 1.0))
matz = MeshPhongMaterial(color = RGBA(0.8, 0.8, 0.8, 1.0))
objx = Rect(Vec(0, 0, 0),Vec(0.50, 0.05, 0.05))
objy = Rect(Vec(0, 0, 0),Vec(0.05, 0.50, 0.05))
objz = Rect(Vec(0, 0, 0),Vec(0.05, 0.05, 0.50))

setobject!(vis2[:world][:x], objx, matx)
setobject!(vis2[:world][:y], objy, maty)
setobject!(vis2[:world][:z], objz, matz)

setobject!(vis2[:body][:x], objx, matx)
setobject!(vis2[:body][:y], objy, maty)
setobject!(vis2[:body][:z], objz, matz)

setobject!(vis2[:vel][:x], objx, matx)
setobject!(vis2[:vel][:y], objy, maty)
setobject!(vis2[:vel][:z], objz, matz)

orientation = :MRP
rt = [0.0, 0.4142138, 0.00]
R = eval(orientation)(rt...)
@show R[:,3]
settransform!(vis2[:body], LinearMap(R))






@variables rt_[1:3]
rt_
R_ = Matrix(MRP(rt_...))
xb_ = R_[:,1]
yb_ = R_[:,2]
zb_ = R_[:,3]

Symbolics.jacobian(xb_, rt_)
Symbolics.jacobian(yb_, rt_)
Symbolics.jacobian(zb_, rt_)







# vis3 = Visualizer()
# open(vis3)
plot_surface!(vis3, model.env)
build_robot!(vis3, model, name=:test)
set_robot!(vis3, model, ref_traj.q[1], name=:test, offset=0.00)
θ = ref_traj.q[1][3]/π

q

Rotations.∇differential(R)
