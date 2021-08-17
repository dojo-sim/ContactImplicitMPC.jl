abstract type Joint end

# masks
x_axis_mask = SMatrix{3,2}([0.0 0.0;
            			    1.0 0.0;
            			    0.0 1.0])

x_axis_input = SMatrix{3,3}([1.0 0.0 0.0;
                			 0.0 0.0 0.0;
                			 0.0 0.0 0.0])

y_axis_mask = SMatrix{3,2}([1.0 0.0;
            			    0.0 0.0;
            			    0.0 1.0])

y_axis_input = SMatrix{3,3}([0.0 0.0 0.0;
            			     0.0 1.0 0.0;
            			     0.0 0.0 0.0])

z_axis_mask = SMatrix{3,2}([1.0 0.0;
            			    0.0 1.0;
            			    0.0 0.0])

z_axis_input = SMatrix{3,3}([0.0 0.0 0.0;
            			     0.0 0.0 0.0;
            			     0.0 0.0 1.0])

# kinematics
function kinematics(r, q)
	# body position
	p = q[1:3]

	# body orientation
	quat = q[4:7]

	k1 = p + quaternion_rotation_matrix(quat) * r

	return k1
end

# rotation axes
function rotation_axes(q_parent, q_child, q_offset)
	# q_parent^-1 q_child q_offset^-1
	(R_multiply(conjugate(q_offset)) * L_multiply(conjugate(q_parent)) * q_child)[2:4]
end

### sphere joint ###
struct SphereJoint{T} <: Joint
    p_parent::SVector{3,T}
    p_child::SVector{3,T}
    p_id::Symbol
    c_id::Symbol
end

dimension(::SphereJoint) = 3

function sphere_constraint(p_parent, p_child, q_parent, q_child)
    kinematics(p_parent, q_parent) - kinematics(p_child, q_child)
end

@variables p_parent[1:3], p_child[1:3], q_parent[1:7], q_child[1:7]

sc = sphere_constraint(p_parent, p_child, q_parent, q_child)
dscp = Symbolics.jacobian(sc, q_parent) * G_func(q_parent)
dscc = Symbolics.jacobian(sc, q_child) * G_func(q_child)

sc_func = eval(Symbolics.build_function(sc, p_parent, p_child, q_parent, q_child)[1])
sc_func! = eval(Symbolics.build_function(sc, p_parent, p_child, q_parent, q_child)[2])

dscp_func = eval(Symbolics.build_function(dscp, p_parent, p_child, q_parent, q_child)[1])
dscp_func! = eval(Symbolics.build_function(dscp, p_parent, p_child, q_parent, q_child)[2])

dscc_func = eval(Symbolics.build_function(dscc, p_parent, p_child, q_parent, q_child)[1])
dscc_func! = eval(Symbolics.build_function(dscc, p_parent, p_child, q_parent, q_child)[2])

_p_parent = rand(3)
_p_child = rand(3)
_q_parent = rand(7)
_q_child = rand(7)
_sc = zeros(3)
_dsc = zeros(3, 6)

@benchmark $_sc .= sc_func($_p_parent, $_p_child, $_q_parent, $_q_child)
@benchmark sc_func!($_sc, $_p_parent, $_p_child, $_q_parent, $_q_child)
@benchmark $_dsc .= dscp_func($_p_parent, $_p_child, $_q_parent, $_q_child)
@benchmark dscp_func!($_dsc, $_p_parent, $_p_child, $_q_parent, $_q_child)
@benchmark $_dsc .= dscc_func($_p_parent, $_p_child, $_q_parent, $_q_child)
@benchmark dscc_func!($_dsc, $_p_parent, $_p_child, $_q_parent, $_q_child)

function constraint(joint::SphereJoint, q_parent, q_child)
    sc_func(joint.p_parent, joint.p_child, q_parent, q_child)
end

function jacobian_parent(joint::SphereJoint, q_parent, q_child)
    dscp_func(joint.p_parent, joint.p_child, q_parent, q_child)
end

function jacobian_child(joint::SphereJoint, q_parent, q_child)
    dscc_func(joint.p_parent, joint.p_child, q_parent, q_child)
end

function constraint!(c, joint::SphereJoint, q_parent, q_child)
    sc_func!(c, joint.p_parent, joint.p_child, q_parent, q_child)
end

function jacobian_parent!(j, joint::SphereJoint, q_parent, q_child)
    dscp_func!(j, joint.p_parent, joint.p_child, q_parent, q_child)
end

function jacobian_child!(j, joint::SphereJoint, q_parent, q_child)
    dscc_func!(j, joint.p_parent, joint.p_child, q_parent, q_child)
end

######

### revolute joint ###
struct RevoluteJoint{T} <: Joint
    p_parent::SVector{3,T}
    p_child::SVector{3,T}
    quat_offset::SVector{4,T}
    mask::SMatrix{3,2,T,6}
    p_id::Symbol
    c_id::Symbol
end

dimension(::RevoluteJoint) = 5

function revolute_constraint(p_parent, p_child, q_parent, q_child, quat_offset, mask)
    [
     kinematics(p_parent, q_parent) - kinematics(p_child, q_child);
     transpose(mask) * rotation_axes(q_parent[4:7], q_child[4:7], quat_offset)
    ]
end

@variables quat_offset[1:4], mask[1:3, 1:2]

rc = revolute_constraint(p_parent, p_child, q_parent, q_child, quat_offset, mask)
drcp = Symbolics.jacobian(rc, q_parent) * G_func(q_parent)
drcc = Symbolics.jacobian(rc, q_child) * G_func(q_child)

rc_func = eval(Symbolics.build_function(rc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[1])
rc_func! = eval(Symbolics.build_function(rc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[2])

drcp_func = eval(Symbolics.build_function(drcp, p_parent, p_child, q_parent, q_child, quat_offset, mask)[1])
drcp_func! = eval(Symbolics.build_function(drcp, p_parent, p_child, q_parent, q_child, quat_offset, mask)[2])

drcc_func = eval(Symbolics.build_function(drcc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[1])
drcc_func! = eval(Symbolics.build_function(drcc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[2])

_quat_offset = rand(4)
_mask = rand(3, 2)
_rc = rand(5)
_drc = rand(5, 6)

@benchmark $_rc .= rc_func($_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark rc_func!($_rc, $_p_parent, $_p_child, $_q_parent, $_q_child,  $_quat_offset, $_mask)
@benchmark $_drc .= drcp_func($_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark drcp_func!($_drc, $_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark $_drc .= drcc_func($_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark drcc_func!($_drc, $_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)

function constraint(joint::RevoluteJoint, q_parent, q_child)
    rc_func(joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_parent(joint::RevoluteJoint, q_parent, q_child)
    drcp_func(joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_child(joint::RevoluteJoint, q_parent, q_child)
    drcc_func(joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function constraint!(c, joint::RevoluteJoint, q_parent, q_child)
    rc_func!(c, joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_parent!(j, joint::RevoluteJoint, q_parent, q_child)
    drcp_func!(j, joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_child!(j, joint::RevoluteJoint, q_parent, q_child)
    drcc_func!(j, joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end
######

### Prismatic Joint ###
struct PrismaticJoint{T} <: Joint
    p_parent::SVector{3,T}
    p_child::SVector{3,T}
    quat_offset::SVector{4,T}
    mask::SMatrix{3,2,T,6}
    p_id::Symbol
    c_id::Symbol
end

dimension(::PrismaticJoint) = 5

function prismatic_constraint(p_parent, p_child, q_parent, q_child, quat_offset, mask)
    [
     transpose(mask) * (quaternion_rotation_matrix(q_parent[4:7]) * kinematics(p_parent, q_parent) - quaternion_rotation_matrix(q_child[4:7]) * kinematics(p_child, q_child))
     rotation_axes(q_parent[4:7], q_child[4:7], quat_offset)
    ]
end

pc = prismatic_constraint(p_parent, p_child, q_parent, q_child, quat_offset, mask)
dpcp = Symbolics.jacobian(pc, q_parent) * G_func(q_parent)
dpcc = Symbolics.jacobian(pc, q_child) * G_func(q_child)

pc_func = eval(Symbolics.build_function(pc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[1])
pc_func! = eval(Symbolics.build_function(pc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[2])

dpcp_func = eval(Symbolics.build_function(dpcp, p_parent, p_child, q_parent, q_child, quat_offset, mask)[1])
dpcp_func! = eval(Symbolics.build_function(dpcp, p_parent, p_child, q_parent, q_child, quat_offset, mask)[2])

dpcc_func = eval(Symbolics.build_function(dpcc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[1])
dpcc_func! = eval(Symbolics.build_function(dpcc, p_parent, p_child, q_parent, q_child, quat_offset, mask)[2])

_pc = rand(5)
_dpc = rand(5, 6)

@benchmark $_pc .= pc_func($_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark pc_func!($_pc, $_p_parent, $_p_child, $_q_parent, $_q_child,  $_quat_offset, $_mask)
@benchmark $_dpc .= dpcp_func($_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark dpcp_func!($_dpc, $_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark $_dpc .= dpcc_func($_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)
@benchmark dpcc_func!($_dpc, $_p_parent, $_p_child, $_q_parent, $_q_child, $_quat_offset, $_mask)

function constraint(joint::PrismaticJoint, q_parent, q_child)
    pc_func(joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_parent(joint::PrismaticJoint, q_parent, q_child)
    dpcp_func(joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_child(joint::PrismaticJoint, q_parent, q_child)
    dpcc_func(joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function constraint!(c, joint::PrismaticJoint, q_parent, q_child)
    pc_func!(c, joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_parent!(joint::PrismaticJoint, q_parent, q_child)
    dpcp_func!(j, joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

function jacobian_child!(joint::PrismaticJoint, q_parent, q_child)
    dpcc_func!(j, joint.p_parent, joint.p_child, q_parent, q_child, joint.quat_offset, joint.mask)
end

######
