################################################################################
# Jacobian of rotated vectors wrt the rotations parameterization.
################################################################################

# Rotation
r = rand(3)
R = MRP(r...)

# Vectors to be rotated
v0 = rand(3)
α = rand()
v1 = α*v0

# Jacobian of rotated vector wrt to r
∇r_v0 = Rotations.∇rotate(R, v0)
∇r_v0_ = ForwardDiff.jacobian(r -> MRP(r...)*v0, r)
@test norm(∇r_v0_ - ∇r_v0) < 1e-8
∇r_v1 = Rotations.∇rotate(R, v1)
@test norm(∇r_v0 - ∇r_v1/α) < 1e-8
