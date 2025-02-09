function _derivative(f, x; h=1e-2)
    return (-f(x + 2h) + 8f(x + h) - 8f(x - h) + f(x - 2h)) / (12 * h)
end

function _second_derivative(f, x; dx=1e-2)
    return (f(x + dx) + f(x - dx) - 2 * f(x)) / dx^2
end

function _sg_filter_second_derivative(f, x; dx=1e-2)
    return (-1 * f(x - 2 * dx) + 16 * f(x - dx) - 30 * f(x) + 16 * f(x + dx) - 1 * f(x + 2 * dx)) / (12 * dx^2)
end

function _first_derivative_smooth(f, x; dx=1e-2)
    return (-2 * f(x - 2 * dx) - f(x + dx) + f(x + dx) + 2 * f(x + 2 * dx)) / (10 * dx)
end
