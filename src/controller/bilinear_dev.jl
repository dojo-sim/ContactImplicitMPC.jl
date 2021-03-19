
name = "quadruped"
model = get_model(name)

nz = num_var(model)
nθ = num_data(model)
z = rand(nz)
θ = rand(nθ)
κ = 1e-5
r0 = rand(nz)
r1 = rand(nz)
r2 = rand(nz)
r3 = rand(nz)
rz0 = spzeros(nz,nz)
rz0 = similar(model.spa.rz_sp, Float64)
rθ0 = zeros(nz,nθ)
lin = LinStep(model, z, θ, κ)

# Test r!
model.res.r(r0, z, θ, κ)
@benchmark r_approx!(lin, r1, z, θ, κ)
@test norm(r0 - r1, Inf) < 1e-8

# Test rz!
model.res.rz(rz0, z, θ)
@benchmark rz_approx!(lin, rz1, z, θ)
rz_FD = rz_approx_FD!(model, lin, rz1, z, θ)
@test norm(rz0 - rz1, Inf) < 1e-8
@test norm(rz1 - rz_FD, Inf) < 1e-8




cg_r_approx!, cg_rz_approx!, cg_rθ_approx! = bilinear_code_gen(model)

lin.z0 = SizedVector{nz}(lin.z0)
lin.θ0 = SizedVector{nθ}(lin.θ0)
lin.r0 = SizedVector{nz}(lin.r0)
lin.rz0
lin.rθ0

rz2 = zeros(nz,nz)
rθ2 = zeros(nz,nθ)
@benchmark cg_r_approx!(r2, z, θ, κ, lin.z0, lin.θ0, lin.r0, lin.rz0, lin.rθ0)
@benchmark cg_rz_approx!(rz2, z, lin.rz0)
@benchmark cg_rθ_approx!(rθ2, lin.rθ0)

cg_r_approx!(r2, z, θ, κ, lin.z0, lin.θ0, lin.r0, lin.rz0, lin.rθ0)
@test norm(r2 - r1, Inf) < 1e-10
cg_rz_approx!(rz2, z, lin.rz0)
@test norm(rz2 - rz1, Inf) < 1e-10

function cg_r_approx2!(r2, z, θ, κ, lin::LinStep)
	cg_r_approx!(r2, z, θ, κ, lin.z0, lin.θ0, lin.r0, lin.rz0, lin.rθ0)
	return nothing
end

@benchmark cg_r_approx2!(r2, z, θ, κ, lin)
