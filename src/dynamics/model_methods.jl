function dynamics(model::QuadrupedModel, dt::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}

	v = (q2[3:end] - q1[3:end]) / dt
	joint_fric = [zeros(2); model.joint_friction * v]

	return (1.0 / dt *
	(
	M_func(model, q0) * (q1 - q0)
	- M_func(model, q1) * (q2 - q1)
	)
	+ transpose(B_func(model, q2)) * u1
	+ transpose(N_func(model, q2)) * γ1
	+ transpose(_P_func(model, q2)) * b1
	- dt * joint_fric
	- dt * _C_func(model, q2, (q2 - q1) / dt))
end


function ∇q0_dynamics(model::QuadrupedModel, dt::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, dt, x, q1, u1, γ1, b1, q2)
	return ForwardDiff.jacobian(dynx, q0)
end

function ∇q1_dynamics(model::QuadrupedModel, dt::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, dt, q0, x, u1, γ1, b1, q2)
	return ForwardDiff.jacobian(dynx, q1)
end

function ∇u1_dynamics(model::QuadrupedModel, dt::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, dt, q0, q1, x, γ1, b1, q2)
	return ForwardDiff.jacobian(dynx, u1)
end

function ∇γ1_dynamics(model::QuadrupedModel, dt::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, dt, q0, q1, u1, x, b1, q2)
	return ForwardDiff.jacobian(dynx, γ1)
end

function ∇b1_dynamics(model::QuadrupedModel, dt::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, dt, q0, q1, u1, γ1, x, q2)
	return ForwardDiff.jacobian(dynx, b1)
end

function ∇q2_dynamics(model::QuadrupedModel, dt::T, q0::Vq0, q1::Vq1,
    u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, dt, q0, q1, u1, γ1, b1, x)
	return ForwardDiff.jacobian(dynx, q2)
end
