
@testset "Quaternion Midpoint" begin
	# Test quaternion sqrt
	q0 = rand(4)
	q0 ./= norm(q0)
	sq0 = sqrt_quat(q0)
	Q0 = Quaternion(q0...)
	sQ0 = sqrt(Q0)
	sq0_ = [sQ0.s, sQ0.v1, sQ0.v2, sQ0.v3]
	@test norm(sq0 - sq0_) < 1e-10

	q0 = [1,0,0,0]
	q0 ./= norm(q0)
	sq0 = sqrt_quat(q0)
	Q0 = Quaternion(q0...)
	sQ0 = sqrt(Q0)
	sq0_ = [sQ0.s, sQ0.v1, sQ0.v2, sQ0.v3]
	@test norm(sq0 - sq0_) < 1e-10

	# Test midpoint
	q0 = rand(4)
	q0 ./= norm(q0)
	ψ = 0.05 * rand(3)
	qψ = φ(ψ)
	sqψ = sqrt_quat(qψ)
	q1 = L_multiply(qψ) * q0
	q2 = L_multiply(sqψ) * L_multiply(sqψ) * q0
	@test norm(q1 - q2) < 1e-10
	@test norm(R_multiply(q0)' * q1 - qψ) < 1e-10

	qmid = L_multiply(sqψ) * q0
	qmid_ = L_multiply(sqrt_quat(R_multiply(q0)' * q1)) * q0
	@test norm(qmid - qmid_) < 1e-10
end
