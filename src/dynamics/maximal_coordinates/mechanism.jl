struct Mechanism{T}
    bodies::Vector{Body}
    sphere_joints::Vector{SphereJoint}
    revolute_joints::Vector{RevoluteJoint}
    prismatic_joints::Vector{PrismaticJoint}
    world_position::SVector{3,T}
    world_orientation::SVector{4,T}
end

function Mechanism()
    Mechanism(Body[],
        SphereJoint[], RevoluteJoint[], PrismaticJoint[],
        SVector{3}([0.0, 0.0, 0.0]), SVector{4}([1.0, 0.0, 0.0, 0.0]))
end

function add!(m::Mechanism, b::Body)
    push!(m.bodies, b)
end

function add!(m::Mechanism, j::SphereJoint)
    push!(m.sphere_joints, j)
end

function add!(m::Mechanism, j::RevoluteJoint)
    push!(m.revolute_joints, j)
end

function add!(m::Mechanism, j::PrismaticJoint)
    push!(m.prismatic_joints, j)
end

#TODO: remove methods
