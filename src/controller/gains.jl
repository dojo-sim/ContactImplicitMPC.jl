function tvlqr(A, B, Q, R)
    T = length(Q)

    P = [zero(A[1]) for t = 1:T]
    K = [zero(B[1]') for t = 1:T-1]
    P[T] = Q[T]

    for t = T-1:-1:1
        K[t] = (R[t] + B[t]' * P[t+1] *  B[t]) \ (B[t]' * P[t+1] * A[t])
        P[t] = (Q[t] + K[t]' * R[t] * K[t]
                + (A[t] - B[t] * K[t])' * P[t+1] * (A[t] - B[t] * K[t]))
    end

    return K, P
end

function reference_gains(model::Model, traj::ContactTraj, obj::TrackingVelocityObjective; N::Int=10, κ=2e-4,
		U_scaling=10, V_scaling=100)
	H = traj.H
	z = [[deepcopy(traj.z) for i = 1:N]...;]
	θ = [[deepcopy(traj.θ) for i = 1:N]...;]

	rz = [LinearizedStep(s, z[t], θ[t], κ).rz for t = 1:N*H]
	rθ = [LinearizedStep(s, z[t], θ[t], κ).rθ for t = 1:N*H]

	idx_z = indices_z(s)
	idx_θ = indices_θ(model, nf=1)

	∂z∂θ = [-rz[t] \ rθ[t] for t = 1:N*H]

	∂q2∂q1 = [zeros(model.nq, model.nq) for t = 1:N*H]
	∂q2∂q2 = [I(model.nq) for t = 1:N*H]
	∂q3∂q1 = [∂z∂θ[t][idx_z.q, idx_θ.q1] for t = 1:N*H]
	∂q3∂q2 = [∂z∂θ[t][idx_z.q, idx_θ.q2] for t = 1:N*H]

	∂q2∂u = [zeros(model.nq, model.nu) for t = 1:N*H]
	∂q3∂u = [∂z∂θ[t][idx_z.q, idx_θ.u] for t = 1:N*H]

	A = [[∂q2∂q1[t] ∂q2∂q2[t]; ∂q3∂q1[t] ∂q3∂q2[t]] for t = 1:N*H]
	B = [[∂q2∂u[t]; ∂q3∂u[t]] for t = 1:N*H]
	Q = [
		[obj.q[1]+V_scaling*obj.v[1] -V_scaling*obj.v[1];
		 -V_scaling*obj.v[1] obj.q[1]+V_scaling*obj.v[1]
		 ] for i = 1:N*H+1]
	R = [U_scaling * obj.u[1] for i = 1:N*H]

	K, P = tvlqr(A, B, Q, R)
	return K[1:H]
end
