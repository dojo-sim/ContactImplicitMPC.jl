function stairs3(x)
    s1 = 0.125
    s2 = 0.375
    s3 = 0.625
    s4 = 0.875
    a0 = 0.00
    a1 = 0.25
    a2 = 0.50
    a3 = 0.75
    y = IfElse.ifelse(x[1] < s1, a0,
            IfElse.ifelse(x[1] < s2, a1,
                IfElse.ifelse(x[1] < s3, a2,
                    IfElse.ifelse(x[1] < s4, a3, a0))))
    return y
end

function d_stairs3(x)
    return 0.0
end

stairs3_2D_lc = Environment{R2, LinearizedCone}(stairs3, d_stairs3)

# X = Vector(0:0.01:1.5)
# plot(X, stairs3steps.(X))

# Kernel
function f(x, center, radius, sharpness)
    return sharpness*(-((x-center)/radius)^2 + 1.0)
end

function smoothed_stairs(x)
    c = [0.0, 0.5, 1.0, 1.5, 2.0]
    a = [0.00, 0.25, 0.50, 0.75, 0.00]
    r = 0.25
    s = 3
    # Kernel vector
    v = [f(x, c[i], r, s) for i=1:5]
    # softmax
    ev = exp.(v)
    w = ev/sum(ev)
    # altitude selection
    return w'*a
end

# X = Vector(0:0.01:3)
# plot(X, smoothed_stairs.(X))
