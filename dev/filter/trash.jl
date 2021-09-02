m = rand(3,3)
m = SMatrix{3,3}(m*m')
# Σ = Diagonal([0.5, 1.0, 2.0])



P1q2 = P1[nq .+ (1:nq), nq .+ (1:nq)]
va = real.(eigvals(sqrt(P1q2)))
ve = real.(eigvecs(sqrt(P1q2)))
ver = copy(ve)
ver[:,1] = -ve[:,1]
det(ve)
det(ver)


vis_3d_gaussian2!(vis, μ, ver, va)

eigvecs(P1)



function ez_rot(θ)
	r = zeros(3,3)
	c = cos(θ)
	s = sin(θ)
	r = [c -s  0;
	     s  c  0;
		 0  0  1]
	return r
end

θ = π/6
R = ez_rot(θ)
ve = eigvecs(R)
va = eigvals(R)
det(ve)
vis_3d_gaussian2!(vis, μ, R, [3,1,0.5])
