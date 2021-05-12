@testset "Environment" begin
	# 2D
	surf = x -> 0.05 * sin(π * x[1])
	surf_grad = x -> [0.05 * π * cos(π * x[1])]
	x0 = ones(10)

	env = ContactControl.environment_2D(surf)
	@test norm(env.surf(x0) - surf(x0)) < 1e-10
	@test norm(env.surf_grad(x0) - surf_grad(x0)) < 1e-10
	@test size(env.surf(x0)) == ()
	@test size(env.surf_grad(x0)) == (1,)

	# 3D
	surf = x -> 0.05 * sin(π * x[1])
	surf_grad = x -> [0.05 * π * cos(π * x[1]); 0.0]
	x0 = ones(10)

	env = ContactControl.environment_3D(surf)
	@test norm(env.surf(x0) - surf(x0)) < 1e-10
	@test norm(env.surf_grad(x0) - surf_grad(x0)) < 1e-10
	@test size(env.surf(x0)) == ()
	@test size(env.surf_grad(x0)) == (2,)
end
