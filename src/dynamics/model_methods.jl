function dynamics(model::QuadrupedModel, h::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}

	v = (q2[3:end] - q1[3:end]) / h[1]
	joint_fric = [zeros(2); model.joint_friction * v]

	return (1.0 / h *
	(
	M_func(model, q0) * (q1 - q0)
	- M_func(model, q1) * (q2 - q1)
	)
	+ transpose(B_func(model, q2)) * u1
	+ transpose(N_func(model, q2)) * γ1
	+ transpose(P_func(model, q2)) * b1
	- h * joint_fric
	- h * C_fast(model, q2, (q2 - q1) / h))
end


function ∇q0_dynamics(model::QuadrupedModel, h::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, h, x, q1, u1, γ1, b1, q2)
	return ForwardDiff.jacobian(dynx, q0)
end

function ∇q1_dynamics(model::QuadrupedModel, h::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, h, q0, x, u1, γ1, b1, q2)
	return ForwardDiff.jacobian(dynx, q1)
end

function ∇u1_dynamics(model::QuadrupedModel, h::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, h, q0, q1, x, γ1, b1, q2)
	return ForwardDiff.jacobian(dynx, u1)
end

function ∇γ1_dynamics(model::QuadrupedModel, h::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, h, q0, q1, u1, x, b1, q2)
	return ForwardDiff.jacobian(dynx, γ1)
end

function ∇b1_dynamics(model::QuadrupedModel, h::T, q0::Vq0, q1::Vq1,
	u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, h, q0, q1, u1, γ1, x, q2)
	return ForwardDiff.jacobian(dynx, b1)
end

function ∇q2_dynamics(model::QuadrupedModel, h::T, q0::Vq0, q1::Vq1,
    u1::Vu1, γ1::Vγ1, b1::Vb1, q2::Vq2) where {T,Vq0,Vq1,Vu1,Vγ1,Vb1,Vq2}
	dynx(x) = dynamics(model, h, q0, q1, u1, γ1, b1, x)
	return ForwardDiff.jacobian(dynx, q2)
end


function residual(model::QuadrupedModel, h::T, z::AbstractVector, θ::AbstractVector, κ) where {T}
	nγ = model.dim.γ
	nb = model.dim.b

	q0, q1, u1 = unpack_θ(model, θ)
	q2, γ1, b1, ψ, η, s1, s2 = unpack_z(model, z)

	ϕ = ϕ_func(model, q2)        # signed-distance function
	# p1 = tangential_contact_pos(model, q1)
	# p2 = tangential_contact_pos(model, q2)
	# vT = (p2 .- p1) ./ model.h
	# vTi = reshape(vT, (Int(nb/nγ),nγ))
	# lin_vT = vcat([[vTi[:,i]; -vTi[:,i]] for i=1:nγ]...)

	lin_vT = P_func(model, q2)*(q2 .- q1) ./ model.h
	# action optimality conditions
	[
	dynamics(model, h, q0, q1, u1, γ1, b1, q2);
	s1 - ϕ;
	γ1 .* s1 .- κ;

	# maximum dissipation optimality conditions
	lin_vT + vcat([ψi.*ones(Int(nb/nγ)) for ψi in ψ]...) - η;
	s2 .- (model.μ * γ1 .- transpose(sum(reshape(b1, (Int(nb/nγ),nγ)), dims=1))[:,1]);
	ψ .* s2 .- κ;
	b1 .* η .- κ
	]
end

function ∇z_residual(model::QuadrupedModel, h::T, z::AbstractVector, θ::AbstractVector, κ) where {T}
	r(x) = residual(model, h, x, θ, κ)
	return ForwardDiff.jacobian(r, z)
end

function ∇θ_residual(model::QuadrupedModel, h::T, z::AbstractVector, θ::AbstractVector, κ) where {T}
	r(x) = residual(model, h, z, x, κ)
	return ForwardDiff.jacobian(r, θ)
end
