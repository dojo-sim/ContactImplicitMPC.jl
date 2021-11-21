@testset "Controller: Linearized Solver" begin
	s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
	model = s.model
	env = s.env

	# Sizes
	nq = model.nq
	nc = model.nc
	nb = nc * friction_dim(env)
	nx = nq
	ny = nb + 2nc
	nz = ContactImplicitMPC.num_var(model, env)
	nθ = ContactImplicitMPC.num_data(model)

	# Indices
	rz_ = ContactImplicitMPC.RZLin(s, rand(nz,nz))
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
	s.res.r!(r0, z0, θ0, [κ])
	s.res.rz!(rz0, z0, θ0)
	s.res.rθ!(rθ0, z0, θ0)

	# Test rz!
	rz1 = ContactImplicitMPC.RZLin(s, rz0)
	ContactImplicitMPC.rz!(rz1, z0, θ0)

	@test norm(rz1.Dx  - rz0[idyn, ix],  Inf) < 1e-10
	@test norm(rz1.Dy1 - rz0[idyn, iy1], Inf) < 1e-10
	@test norm(rz1.Rx  - rz0[irst, ix],  Inf) < 1e-10
	@test norm(rz1.Ry1 - rz0[irst, iy1],  Inf) < 1e-10
	@test norm(rz1.Ry2 - diag(rz0[irst, iy2]),  Inf) < 1e-10
	@test norm(rz1.y2  - diag(rz0[ibil, iy1]),  Inf) < 1e-10
	@test norm(rz1.y1  - diag(rz0[ibil, iy2]),  Inf) < 1e-10

	# Test r!
	r1 = ContactImplicitMPC.RLin(s, z0, θ0, r0, rz0, rθ0)
	ContactImplicitMPC.rlin!(r1, z0, θ0, [κ])

	@test norm(r0[r1.idyn] - r1.rdyn, Inf) < 1e-10
	@test norm(r0[r1.irst] - r1.rrst, Inf) < 1e-10
	@test norm(r0[r1.ibil] - r1.rbil, Inf) < 1e-10

	# Test linear_solve!
	Δ = rand(nz)
	ContactImplicitMPC.linear_solve!(Δ, rz1, r1)
	@test norm(Δ - rz0 \ r0, Inf) < 1e-10

	# Test linear_solve! for matrices
	s.res.rθ!(rθ0, z0, θ0)
	@test norm(rθ0[ibil, :], Inf) < 1e-10

	rz1 = ContactImplicitMPC.RZLin(s, rz0)
	rθ1 = ContactImplicitMPC.RθLin(s, rθ0)
	δz1 = rand(nz,nθ)
	ContactImplicitMPC.linear_solve!(δz1, rz1, rθ1)
	@test norm(δz1 - (rz0 \ rθ0), Inf) < 1e-10
end

@testset "Controller: Update Linearized Residuals" begin
	s = get_simulation("hopper_2D", "flat_2D_lc", "flat")
	model = s.model
	env = s.env
	nz = ContactImplicitMPC.num_var(model, env)
	nθ = ContactImplicitMPC.num_data(model)

	z0 = rand(nz)
	z1 = rand(nz)
	θ0 = rand(nθ)
	θ1 = rand(nθ)
	rz0 = zeros(nz,nz)
	rz1 = zeros(nz,nz)
	s.res.rz!(rz0, z0, θ0)
	s.res.rz!(rz1, z1, θ1)
	rz = ContactImplicitMPC.RZLin(s, rz0)
	ContactImplicitMPC.update!(rz, rz1)
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
	r = ContactImplicitMPC.RLin(s, z0, θ0, r0, rz0, rθ0)
	ContactImplicitMPC.update!(r, z1, θ1, r1, rz1, rθ1)

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
