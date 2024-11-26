"""
    terms(f::Function,Λ::Float64,n::Int64)

Returns the terms for the separable representation of the function `f` given by
f(p,p') = ∑ᵢ gᵢ(p)gᵢ(p') where gᵢ are the terms returned by this function.

# Arguments
- `f::Function`: The function to be represented.
- `Λ::Float64`: Momentum cutoff.
- `n::Int64`: Number of terms in approximation to be returned.
"""
function _terms(f::Function, points, n)
    func_list = []
    cur_F(x, y) = f(x, y)
    push!(func_list, cur_F)
    for j = 1:(n-1)
        cur_F = next_func(cur_F, points[j])
        push!(func_list, cur_F)
    end
    return red_arg_list(func_list, points)
end

function terms(f::Function, n::Int64, Λ::Float64)
    points = collect(range(0, Λ, length=n))
    return _terms(f, points, n)
end

function terms(f::Function, points)
    n = length(points)
    return _terms(f, points, n)
end

function terms(f::Function, n::Int64)
    return nothing
end

## Helper functions
function next_func(func::Function, p0)
    g(x, y) = func(x, y) - func(x, p0) * func(p0, y) / func(p0, p0)
    return g
end

function red_arg(f::Function, point::Float64)
    return x -> f(x, point) / sqrt(f(point, point))
end

function red_arg_list(f_list, points)
    return [red_arg(f_list[i], points[i]) for i = 1:length(f_list)]
end

export terms
