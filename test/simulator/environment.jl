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


# 2D
α = 0.05
surf = x -> α * sin(π * x[1])
surf_grad = x -> [α * π * cos(π * x[1]); ]
x0 = ones(10)
env = ContactControl.environment_2D(surf)
@test norm(env.surf(x0) - surf(x0)) < 1e-10
@test norm(env.surf_grad(x0) - surf_grad(x0)) < 1e-10
@test size(env.surf(x0)) == ()
@test size(env.surf_grad(x0)) == (1,)
@test norm(env.surf(x0) - 0.0) < 1e-10
@test norm(env.surf_grad(x0) + [α*π]) < 1e-10


verify_2D_surface(env, x0; x_range = range(0.0, stop = 5.0, length = 1000))


# 3D
α = 1.0
surf = x -> α * (sin(π * x[1]) + sin(π * x[2]))
surf_grad = x -> α * [ π * cos(π * x[1]),  0.0*π * cos(π * x[2])]
x0 = ones(10)

env = ContactControl.environment_3D(surf)
@test norm(env.surf(x0) - surf(x0)) < 1e-10
@test norm(env.surf_grad(x0) - surf_grad(x0)) < 1e-10
@test size(env.surf(x0)) == ()
@test size(env.surf_grad(x0)) == (2,)
@test norm(env.surf(x0) - 0.0) < 1e-10
@test norm(env.surf_grad(x0) + [α*π, α*π]) < 1e-10

verify_3D_surface(env, [x0]; x_range = range(0.0, stop = 5.0, length = 1000))

env = environment_3D_flat()
rotation(env, rand(10))


const ContactControl = Main
