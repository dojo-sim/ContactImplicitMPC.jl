include(joinpath(pwd(), "src/dynamics/maximal_coordinates/utils.jl"))
include(joinpath(pwd(), "src/dynamics/maximal_coordinates/body.jl"))
include(joinpath(pwd(), "src/dynamics/maximal_coordinates/joint.jl"))
include(joinpath(pwd(), "src/dynamics/maximal_coordinates/mechanism.jl"))

mechanism = Mechanism()

body1 = Rod(1.0, Diagonal(@SVector ones(3)), id = :body1)
body2 = Rod(1.0, Diagonal(@SVector ones(3)), id = :body2)

joint1 = SphereJoint(
    SVector{3}([0.0, 0.0, 0.0]),
    SVector{3}([0.0, 0.0, 0.0]),
    body1.id,
    body2.id)

push!(mechanism.bodies, body1)
push!(mechanism.bodies, body2)
push!(mechanism.joints, joint1)

vis = Visualizer()
render(vis)
visualize!(vis, mechanism, [])
