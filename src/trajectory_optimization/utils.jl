"""
    linear interpolation between two vectors over horizon T
"""
function linear_interpolation(x0, xf, T)
    n = length(x0)
    X = [copy(Array(x0)) for t = 1:T]
    for t = 1:T
        for i = 1:n
            X[t][i] = (xf[i] - x0[i]) / (T - 1) * (t - 1) + x0[i]
        end
    end
    return X
end
