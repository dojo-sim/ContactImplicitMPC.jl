quadratic_bowl_3D_lc = environment_3D(x -> transpose(x[1:2]) * x[1:2], cone = LinearizedCone)
quadratic_bowl_3D_nc = environment_3D(x -> transpose(x[1:2]) * x[1:2], cone = NonlinearCone)

circular_bowl_3D_nc = environment_3D(x -> -sqrt(2.5^2.0 - x[1]^2.0 - x[2]^2.0) + 2.5, cone = NonlinearCone)

# r = 2.0
# t = range(-r, stop = r, length = 1000)
# z = -sqrt.(r^2.0 .- t.^2.0 .- 0^2.0) .+ r
#
# plot(t, z, aspect_ratio = :equal)
# env = deepcopy(circular_bowl_3D_nc)
# env = deepcopy(flat_3D_nc)
# q = [-0.5; -0.5; env.surf([-0.5; -0.5])]
# n = [env.surf_grad(q[1:2]); -1.0]
# ns = n ./ sqrt(transpose(n) * n) #* model.r
