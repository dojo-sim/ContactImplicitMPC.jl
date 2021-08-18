include(joinpath(pwd(), "src/dynamics/maximal_coordinates/utils.jl"))
include(joinpath(pwd(), "src/dynamics/maximal_coordinates/body.jl"))
include(joinpath(pwd(), "src/dynamics/maximal_coordinates/joint.jl"))
include(joinpath(pwd(), "src/dynamics/maximal_coordinates/mechanism.jl"))
include(joinpath(pwd(), "src/dynamics/maximal_coordinates/visualize.jl"))

mechanism = Mechanism()

body1 = Rod(1.0, Diagonal(@SVector ones(3)), id = :body1)
body2 = Rod(1.0, Diagonal(@SVector ones(3)), id = :body2)

joint1 = SphereJoint(
    SVector{3}([0.0, 0.0, 0.0]),
    SVector{3}([0.0, 0.0, 0.0]),
    body1.id,
    body2.id)

add!(mechanism, body1)
add!(mechanism, body2)
add!(mechanism, joint1)

position1 = SVector{3}([0.0, 0.0, -0.5 * body1.geometry[:length]])
orientation1 = static_quaternion(one(UnitQuaternion))
q1 = [position1; orientation1]

position2 = SVector{3}([0.0, 0.0, -body1.geometry[:length] - 0.5 * body2.geometry[:length]])
orientation2 = static_quaternion(one(UnitQuaternion))
q2 = [position2; orientation2]

configurations = [q1, q2]

vis = Visualizer()
render(vis)
create_mechanism!(vis, mechanism)
configuration!(vis, mechanism, configurations)

mechanism
