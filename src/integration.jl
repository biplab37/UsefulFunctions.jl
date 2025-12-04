## This file contains some of the integration routines that I use frequently.

function trapizoid(f::Function, a::Float64, b::Float64, n::Int64)
    """
    trapizoid(f::Function, a::Float64, b::Float64, n::Int64)

    Returns the integral of `f` from `a` to `b` using the trapezoidal rule with `n` steps.
    """
    h = (b-a)/n
    sum = 0.5*(f(a) + f(b))
    for i in 1:n-1
        sum += f(a + i*h)
    end
    return h*sum
end

function simpson(f::Function, a::Float64, b::Float64, n::Int64)
    """
    simpson(f::Function, a::Float64, b::Float64, n::Int64)

    Returns the integral of `f` from `a` to `b` using Simpson's rule with `n` steps.
    """
    h = (b-a)/n
    sum = f(a) + f(b)
    for i in 1:n-1
        if i%2 == 0
            sum += 2*f(a + i*h)
        else
            sum += 4*f(a + i*h)
        end
    end
    return h*sum/3
end

export trapizoid, simpson