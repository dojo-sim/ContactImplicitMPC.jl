
const ContactControl = Main
s = get_simulation("quadruped", "flat_2D_lc", "flat")
model = s.model
env = s.env

# Sizes
nq = model.dim.q
nc = model.dim.c
nb = nc * friction_dim(env)
nx = nq
ny = nb + 2nc
nz = ContactControl.num_var(model, env)
nθ = ContactControl.num_data(model)

# Indices
rz_ = ContactControl.RZLin(s, rand(nz,nz))
ix = rz_.ix
iy1 = rz_.iy1
iy2 = rz_.iy2
idyn = rz_.idyn
irst = rz_.irst
ibil = rz_.ibil

z0 = rand(nz)
θ0 = rand(nθ)
r0 = rand(nz)
κ = 1e-4
rz0 = zeros(nz,nz)
rθ0 = zeros(nz,nθ)
s.res.r!(r0, z0, θ0, κ)
s.res.rz!(rz0, z0, θ0)
s.res.rθ!(rθ0, z0, θ0)
vix = Vector(rz_.ix)
viy1 = Vector(rz_.iy1)
viy2 = Vector(rz_.iy2)
vidyn = Vector(rz_.idyn)
virst = Vector(rz_.irst)
vibil = Vector(rz_.ibil)
plot(Gray.(1e10*abs.(rz0[[vidyn; virst; vibil;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz0[[vidyn;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz0[[virst;], [vix; viy1; viy2]])))
plot(Gray.(1e10*abs.(rz0[[vibil;], [vix; viy1; viy2]])))


# Test rz!
rz1 = ContactControl.RZLin(s, rz0)
ContactControl.rz!(rz1, z0)

@test norm(rz1.Dx  - rz0[idyn, ix],  Inf) < 1e-10
@test norm(rz1.Dy1 - rz0[idyn, iy1], Inf) < 1e-10
@test norm(rz1.Rx  - rz0[irst, ix],  Inf) < 1e-10
@test norm(rz1.Ry1 - rz0[irst, iy1],  Inf) < 1e-10
@test norm(rz1.Ry2 - diag(rz0[irst, iy2]),  Inf) < 1e-10
@test norm(rz1.y2  - diag(rz0[ibil, iy1]),  Inf) < 1e-10
@test norm(rz1.y1  - diag(rz0[ibil, iy2]),  Inf) < 1e-10

# Test r!
r1 = ContactControl.RLin(s, z0, θ0, r0, rz0, rθ0)
ContactControl.r!(r1, z0, θ0, κ)
rz1.S.D
plot(Gray.(1e10*abs.(Matrix(rz1.S.D))))

@test norm(r0[r1.idyn] - r1.rdyn, Inf) < 1e-10
@test norm(r0[r1.irst] - r1.rrst, Inf) < 1e-10
@test norm(r0[r1.ibil] - r1.rbil, Inf) < 1e-10

# Test linear_solve!
Δ = rand(nz)
ContactControl.linear_solve!(Δ, rz1, r1)
@test norm(Δ - rz0 \ r0, Inf) < 1e-10

# @benchmark ContactControl.rz!(rz1, z)
# @benchmark ContactControl.r!(r1, z, θ, κ)
# @benchmark ContactControl.linear_solve!(Δ, rz1, r1)

# Test linear_solve! for matrices
s.res.rθ!(rθ0, z0, θ0)
@test norm(rθ0[ibil, :], Inf) < 1e-10
plot(Gray.(1e10*abs.(rθ0[idyn,:])))
plot(Gray.(1e10*abs.(rθ0[irst,:])))
plot(Gray.(1e10*abs.(rθ0[ibil,:])))

rz1 = ContactControl.RZLin(s, rz0)
rθ1 = ContactControl.RθLin(s, rθ0)
δz1 = rand(nz,nθ)
ContactControl.linear_solve!(δz1, rz1, rθ1)
@test norm(δz1 - (rz0 \ rθ0), Inf) < 1e-10
# @benchmark ContactControl.linear_solve!(δz1, rz1, rθ1)

using JLD2
path = joinpath(module_dir(), "dev", "rz.jld2")
ordered_rz = rz0[[vidyn; virst; vibil;], [vix; viy1; viy2]]

@save path ordered_rz
