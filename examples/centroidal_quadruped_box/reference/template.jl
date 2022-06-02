function ellipse_traj(x_start, x_goal, z, T)
	dist = x_goal - x_start
	a = 0.5 * dist
	b = z
	z̄ = 0.0
	# x = range(x_start, stop = x_goal, length = T)
	x = circular_projection_range(x_start, stop = x_goal, length = T)
	z = sqrt.(max.(0.0, (b^2) * (1.0 .- ((x .- (x_start + a)).^2.0) / (a^2.0))))
	return x, z
end

function circular_projection_range(start; stop=1.0, length=10)
	dist = stop - start
	θr = range(π, stop=0, length=length)
	r = start .+ dist * ((1 .+ cos.(θr))./2)
	return r
end