function _derivative(f, x; h=1e-2)
    return (-f(x + 2h) + 8f(x + h) - 8f(x - h) + f(x - 2h)) / (12 * h)
end

export _derivative