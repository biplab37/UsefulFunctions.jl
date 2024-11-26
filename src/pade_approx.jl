Base.@kwdef struct Pade
    m::Int
    n::Int
    P::Vector{Float64}
    Q::Vector{Float64}
    f_approx::Function
end

## pade approximation
"""
    pade_approximation(x_values, f_values, m, n)

Returns pade approximation given two vectors, the x_values and the function evaluated at those points f_values
"""
function pade_approximation(x_values, f_values, m, n)
    @assert length(x_values) == length(f_values) >= m + n + 1
    A = zeros(length(x_values), m + n + 1)
    for (i, k) in enumerate(x_values)
        for j in 1:m+1
            A[i, j] = k^(j - 1)
        end
        for j in 1:n
            A[i, m+j+1] = -f_values[i] * k^j
        end
    end
    coeffs = A \ f_values
    P = coeffs[1:m+1]
    Q = coeffs[m+2:end]
    f(x) = sum([P[i] * x^(i - 1) for i in 1:m+1]) / (1 + sum([Q[i] * x^i for i in 1:n]))

    return Pade(m=m, n=n, P=P, Q=Q, f_approx=f)
end

function pade_approximation(f::Function, m, n, interval)
    a, b = interval
    num = m + n + 1
    x_values = @. (a + b) / 2 + (b - a) * cos((π * (2 * (1:n) - 1)) / (2num)) / 2
    f_values = f.(x_values)
    return pade_approximation(x_values, f_values, m, n)
end

export pade_approximation
