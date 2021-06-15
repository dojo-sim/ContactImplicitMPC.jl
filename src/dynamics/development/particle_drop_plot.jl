using Plots

a = -1.0

t = range(0, stop = 4.0, length = 1000)

a_traj = Float64[]
v_traj = Float64[]
p_traj = Float64[]

flag = false
for tt in t
	p = 2.0 + 0.5 * a * tt^2.0
	p = p <= 0.0 ? 0.0 : p
	push!(p_traj, p)
	v = p <= 0.0 ? 0.0 : a * tt
	push!(v_traj, v)

	# if !flag && p <= 0.0
	# 	flag = true
	# 	push!(a_traj, 10.0)
	# else
	push!(a_traj, a)
	# end
end

plot!(t, p_traj, color = :cyan, width = 4.0, label = "pos.")
plot!(t, v_traj, color = :orange, width = 4.0, label = "vel.")
plot(t, a_traj, color = :magenta, width = 4.0, label = "acc.", xlabel = "time")

using PGFPlots
const PGF = PGFPlots

aa = PGF.Plots.Linear(t, a_traj, mark="none",
	style="color=magenta, line width=2pt", legendentry = "acc.")
aa2 = PGF.Plots.Linear([2.0, 2.0], [a; 10.0], mark="none",
	style="color=magenta, line width=2pt")
vv = PGF.Plots.Linear(t, v_traj, mark="none",
	style="color=orange, line width=2pt", legendentry = "vel.")
pp = PGF.Plots.Linear(t, p_traj, mark="none",
	style="color=cyan, line width=2pt", legendentry = "pos.")

plt = Axis([aa2;aa;vv;pp],
    axisEqualImage=false,
    hideAxis=false,
	xlabel="time",
	ymin = -2.0,
	ymax = 2.0,
	xmin = t[1],
	xmax = t[end],
	legendStyle="{at={(0.95,0.95)},anchor=north east}")

# Save to tikz format
dir = joinpath(pwd(),"src/dynamics/development")
PGF.save(joinpath(dir,"drop_plot.tikz"), plt, include_preamble=false)
