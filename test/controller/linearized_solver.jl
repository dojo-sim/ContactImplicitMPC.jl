# @testset "Controller: Linearized Solver" begin
# 	s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
# 	model = s.model
# 	env = s.env
#
# 	# Sizes
# 	nq = model.dim.q
# 	nc = model.dim.c
# 	nb = nc * friction_dim(env)
# 	nx = nq
# 	ny = nb + 2nc
# 	nz = ContactControl.num_var(model, env)
# 	nθ = ContactControl.num_data(model)
#
# 	# Indices
# 	rz_ = ContactControl.RZLin(s, rand(nz,nz))
# 	ix = rz_.ix
# 	iy1 = rz_.iy1
# 	iy2 = rz_.iy2
# 	idyn = rz_.idyn
# 	irst = rz_.irst
# 	ibil = rz_.ibil
#
# 	z0 = rand(nz)
# 	θ0 = rand(nθ)
# 	r0 = rand(nz)
# 	rz0 = zeros(nz,nz)
# 	rθ0 = zeros(nz,nθ)
# 	s.res.rz!(rz0, z0, θ0)
# 	s.res.rθ!(rθ0, z0, θ0)
# 	# vix = Vector(rz_.ix)
# 	# viy1 = Vector(rz_.iy1)
# 	# viy2 = Vector(rz_.iy2)
# 	# vidyn = Vector(rz_.idyn)
# 	# virst = Vector(rz_.irst)
# 	# vibil = Vector(rz_.ibil)
# 	# plot(Gray.(1e10*abs.(rz0[[vidyn; virst; vibil;], [vix; viy1; viy2]])))
# 	# plot(Gray.(1e10*abs.(rz0[[vidyn;], [vix; viy1; viy2]])))
# 	# plot(Gray.(1e10*abs.(rz0[[virst;], [vix; viy1; viy2]])))
# 	# plot(Gray.(1e10*abs.(rz0[[vibil;], [vix; viy1; viy2]])))
#
# 	z = rand(nz)
# 	θ = rand(nθ)
# 	κ = 1e-4
# 	κv = [κ]
# 	rz1 = ContactControl.RZLin(s, rz0)
# 	r1 = ContactControl.RLin(s, z0, θ0, r0, rz0, rθ0)
#
# 	# Test rz!
# 	ContactControl.rz!(rz1, z)
# 	rz2 = rand(nz, nz)
# 	# s.linearized.rz!(rz2, z, rz0)
# 	#
# 	# @test norm(rz1.Dx  - rz2[idyn, ix],  Inf) < 1e-10
# 	# @test norm(rz1.Dy1 - rz2[idyn, iy1], Inf) < 1e-10
# 	# @test norm(rz1.Rx  - rz2[irst, ix],  Inf) < 1e-10
# 	# @test norm(rz1.Ry1 - rz2[irst, iy1],  Inf) < 1e-10
# 	# @test norm(rz1.Ry2 - diag(rz2[irst, iy2]),  Inf) < 1e-10
# 	# @test norm(rz1.y2  - diag(rz2[ibil, iy1]),  Inf) < 1e-10
# 	# @test norm(rz1.y1  - diag(rz2[ibil, iy2]),  Inf) < 1e-10
#
# 	# Test r!
# 	r1 = ContactControl.RLin(s, z0, θ0, r0, rz0, rθ0)
# 	ContactControl.r!(r1, z, θ, κ)
# 	# r2  = rand(nz)
# 	# model.linearized.r!(r2, z, θ, κv, z0, θ0, r0, rz0, rθ0)
# 	# @test norm(r2[r1.idyn] - r1.rdyn, Inf) < 1e-10
# 	# @test norm(r2[r1.irst] - r1.rrst, Inf) < 1e-10
# 	# @test norm(r2[r1.ibil] - r1.rbil, Inf) < 1e-10
#
# 	# Test linear_solve!
# 	# Δ = rand(nz)
# 	# ContactControl.linear_solve!(Δ, rz1, r1)
# 	# @test norm(Δ - rz2 \ r2, Inf) < 1e-10
# 	# #
# 	# # @benchmark ContactControl.rz!(rz1, z)
# 	# # @benchmark ContactControl.r!(r1, z, θ, κ)
# 	# # @benchmark ContactControl.linear_solve!(Δ, rz1, r1)
# 	#
# 	# # Test linear_solve! for matrices
# 	# model.res.rθ!(rθ0, z, θ)
# 	# @test norm(rθ0[ibil, :], Inf) < 1e-10
# 	# # plot(Gray.(1e10*abs.(rθ_[idyn,:])))
# 	# # plot(Gray.(1e10*abs.(rθ_[irst,:])))
# 	# # plot(Gray.(1e10*abs.(rθ_[ibil,:])))
# 	#
# 	# rz1 = ContactControl.RZLin(model, rz0)
# 	# rθ1 = ContactControl.RθLin(model, rθ0)
# 	# δz1 = rand(nz,nθ)
# 	# ContactControl.linear_solve!(δz1, rz1, rθ1)
# 	# @test norm(δz1 - (rz0 \ rθ0), Inf) < 1e-10
# 	# # @benchmark ContactControl.linear_solve!(δz1, rz1, rθ1)
# end

@testset "Controller: Update Linearized Residuals" begin
	s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
	model = s.model
	env = s.env
	nz = ContactControl.num_var(model, env)
	nθ = ContactControl.num_data(model)

	z0 = rand(nz)
	z1 = rand(nz)
	θ0 = rand(nθ)
	θ1 = rand(nθ)
	rz0 = zeros(nz,nz)
	rz1 = zeros(nz,nz)
	s.res.rz!(rz0, z0, θ0)
	s.res.rz!(rz1, z1, θ1)
	rz = ContactControl.RZLin(s, rz0)
	ContactControl.update!(rz, rz1)
	@test norm(rz.Dx  - rz1[rz.idyn, rz.ix])  < 1e-10
	@test norm(rz.Dy1 - rz1[rz.idyn, rz.iy1]) < 1e-10
	@test norm(rz.Rx  - rz1[rz.irst, rz.ix])  < 1e-10
	@test norm(rz.Ry1 - rz1[rz.irst, rz.iy1]) < 1e-10
	@test norm(rz.Ry2 - diag(rz1[rz.irst, rz.iy2])) < 1e-10
	@test norm(rz.y1  - diag(rz1[rz.ibil, rz.iy2])) < 1e-10
	@test norm(rz.y2  - diag(rz1[rz.ibil, rz.iy1])) < 1e-10

	z0 = rand(nz)
	θ0 = rand(nθ)
	r0 = rand(nz)
	rz0 = zeros(nz,nz)
	rθ0 = zeros(nz,nθ)
	z1 = rand(nz)
	θ1 = rand(nθ)
	r1 = rand(nz)
	rz1 = zeros(nz,nz)
	rθ1 = zeros(nz,nθ)
	s.res.rz!(rz0, z0, θ0)
	s.res.rθ!(rθ0, z0, θ0)
	s.res.rz!(rz1, z1, θ1)
	s.res.rθ!(rθ1, z1, θ1)
	r = ContactControl.RLin(s, z0, θ0, r0, rz0, rθ0)
	ContactControl.update!(r, z1, θ1, r1, rz1, rθ1)

	@test norm(r.rdyn0 - r1[r.idyn])  < 1e-10
	@test norm(r.rrst0 - r1[r.irst]) < 1e-10
	@test norm(r.rbil0 - r1[r.ibil]) < 1e-10

	@test norm(r.Dx  - rz1[rz.idyn, rz.ix])  < 1e-10
	@test norm(r.Dy1 - rz1[rz.idyn, rz.iy1]) < 1e-10
	@test norm(r.Rx  - rz1[rz.irst, rz.ix])  < 1e-10
	@test norm(r.Ry1 - rz1[rz.irst, rz.iy1]) < 1e-10
	@test norm(r.Ry2 - diag(rz1[rz.irst, rz.iy2])) < 1e-10

	@test norm(r.rθdyn - rθ1[r.idyn,:]) < 1e-10
	@test norm(r.rθrst - rθ1[r.irst,:]) < 1e-10
	@test norm(r.rθbil - rθ1[r.ibil,:]) < 1e-10

	@test norm(r.x0  - z1[r.ix])  < 1e-10
	@test norm(r.y10 - z1[r.iy1]) < 1e-10
	@test norm(r.y20 - z1[r.iy2]) < 1e-10
	@test norm(r.θ0 -  θ1) < 1e-10
end



s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
s = get_simulation("quadruped", "flat_2D_lc", "flat")
s = get_simulation("flamingo", "flat_2D_lc", "flat")
s = get_simulation("pushbot", "flat_2D_lc", "flat")
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
rz0 = zeros(nz,nz)
rθ0 = zeros(nz,nθ)
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

plot(Gray.(1e0*abs.(rθ0[[vidyn; virst; vibil;], :])))
plot(Gray.(1e10*abs.(rθ0[[vidyn; virst; vibil;], :])))
plot(Gray.(1e10*abs.(rθ0[[vidyn;], :])))
plot(Gray.(1e10*abs.(rθ0[[virst;], :])))
plot(Gray.(1e10*abs.(rθ0[[vibil;], :])))
nθ
nq + nc + nb

nc
nb
nz
nq + 4nc + 2nb


z = rand(nz)
θ = rand(nθ)
κ = 1e-4
κv = [κ]
rz1 = ContactControl.RZLin(s, rz0)
r1 = ContactControl.RLin(s, z0, θ0, r0, rz0, rθ0)

# Test rz!
ContactControl.rz!(rz1, z)
rz2 = rand(nz, nz)
s.linearized.rz!(rz2, z, rz0)

@test norm(rz1.Dx  - rz2[idyn, ix],  Inf) < 1e-10
@test norm(rz1.Dy1 - rz2[idyn, iy1], Inf) < 1e-10
@test norm(rz1.Rx  - rz2[irst, ix],  Inf) < 1e-10
@test norm(rz1.Ry1 - rz2[irst, iy1],  Inf) < 1e-10
@test norm(rz1.Ry2 - diag(rz2[irst, iy2]),  Inf) < 1e-10
@test norm(rz1.y2  - diag(rz2[ibil, iy1]),  Inf) < 1e-10
@test norm(rz1.y1  - diag(rz2[ibil, iy2]),  Inf) < 1e-10

# Test r!
r1 = ContactControl.RLin(s, z0, θ0, r0, rz0, rθ0)
ContactControl.r!(r1, z, θ, κ)
# r2  = rand(nz)
# model.linearized.r!(r2, z, θ, κv, z0, θ0, r0, rz0, rθ0)
# @test norm(r2[r1.idyn] - r1.rdyn, Inf) < 1e-10
# @test norm(r2[r1.irst] - r1.rrst, Inf) < 1e-10
# @test norm(r2[r1.ibil] - r1.rbil, Inf) < 1e-10

# Test linear_solve!
# Δ = rand(nz)
# ContactControl.linear_solve!(Δ, rz1, r1)
# @test norm(Δ - rz2 \ r2, Inf) < 1e-10
# #
# # @benchmark ContactControl.rz!(rz1, z)
# # @benchmark ContactControl.r!(r1, z, θ, κ)
# # @benchmark ContactControl.linear_solve!(Δ, rz1, r1)
#
# # Test linear_solve! for matrices
# model.res.rθ!(rθ0, z, θ)
# @test norm(rθ0[ibil, :], Inf) < 1e-10
# # plot(Gray.(1e10*abs.(rθ_[idyn,:])))
# # plot(Gray.(1e10*abs.(rθ_[irst,:])))
# # plot(Gray.(1e10*abs.(rθ_[ibil,:])))
#
# rz1 = ContactControl.RZLin(model, rz0)
# rθ1 = ContactControl.RθLin(model, rθ0)
# δz1 = rand(nz,nθ)
# ContactControl.linear_solve!(δz1, rz1, rθ1)
# @test norm(δz1 - (rz0 \ rθ0), Inf) < 1e-10
# # @benchmark ContactControl.linear_solve!(δz1, rz1, rθ1)
