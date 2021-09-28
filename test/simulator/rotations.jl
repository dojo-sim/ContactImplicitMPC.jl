@testset "Simulator: Rotations" begin
	# 3D
	a = cos(0.25 * π)
	env = ContactImplicitMPC.Environment{R3,FrictionCone}(x -> 1.0 * x[1], x -> [1.0; 0.0])

	R = rotation(env, [0.0; 0.0])
	λ = [0.0; 0.0; 1.0]
	@test norm(R' * λ - [-a; 0.0; a]) < 1.0e-8

	λ = [1.0; 0.0; 0.0]
	@test norm(R' * λ - [a; 0.0; a]) < 1.0e-8

	λ = [0.0; 1.0; 0.0]
	@test norm(R' * λ - [0.0; 1.0; 0.0]) < 1.0e-8

	env = ContactImplicitMPC.Environment{R3,FrictionCone}(x -> x' * x, x -> 2.0 * x)

	R = rotation(env, [0.0; 0.0])
	λ = [0.0; 0.0; 1.0]
	@test norm(R' * λ - [0.0; 0.0; 1.0]) < 1.0e-8

	λ = [1.0; 0.0; 0.0]
	@test norm(R' * λ - [1.0; 0.0; 0.0]) < 1.0e-8

	λ = [0.0; 1.0; 0.0]
	@test norm(R' * λ - [0.0; 1.0; 0.0]) < 1.0e-8

	a = 1 / 3
	b = 2 / 3
	R = rotation(env, [1.0; 1.0])
	λ = [0.0; 0.0; 1.0]
	@test norm(R' * λ - [-b; -b; a]) < 1.0e-8

	R' * λ
	λ = [1.0; 0.0; 0.0]
	@test norm(R' * λ - [b; -a; b]) < 1.0e-8

	λ = [0.0; 1.0; 0.0]
	@test norm(R' * λ - [-a; b; b]) < 1.0e-8

	# 2D
	a = cos(0.25 * π)
	x1 = [1.0; 0.0]
	x2 = [0.0; 1.0]
	env = ContactImplicitMPC.Environment{R2,FrictionCone}(x -> 1.0 * x[1], x -> 1.0)

	R = rotation(env, [0.0])
	λ = [0.0; 1.0]
	@test norm(R' * λ - [-a; a]) < 1.0e-8

	λ = [1.0; 0.0]
	@test norm(R' * λ - [a; a]) < 1.0e-8
end
