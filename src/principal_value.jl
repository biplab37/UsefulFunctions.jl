## Contains the function for calculating Principal Value integration.

"""
    PVintegral(f::Function, a, b, c, integrate::Function)

Calculates the Cauchy principal value integral
∫ₐᵇf(x)/(x-c)dx
First it turns into a symmetric integral from -1 to 1, then substracts the singulatity.
It uses the function `integrate` to calculate all the finite integrals.
"""
function PVintegral(f::Function, a, b, c, integrate::Function)
    if a < c < b
        return _principal_asymmetric(f, a, b, c, integrate)
    else
        return integrate(x -> f(x) / (x - c), a, b)
    end
end

function _symmetrize(f::Function, a, c)
    return x -> f((c - a)x + c)
end

function _sub_singularity(f::Function)
    fsing(t) = (t == 0) ? _derivative(f, 0.0) : (f(t) - f(0)) / t
    return fsing
end

function _principal_asymmetric(f, a, b, c, integrate)
    fsym = _symmetrize(f, a, c)
    fsing = _sub_singularity(fsym)
    return integrate(x -> (f(x) - f(c)) / (x - c), 2c - a, b) + f(c) * log((b - c) / (c - a)) + integrate(fsing, -1.0, 1.0)
end

export PVintegral
