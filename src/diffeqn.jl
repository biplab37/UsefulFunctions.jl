@doc raw"""Solves (system of) Ordinary differential equations using the Runge-Kutta method.
    function rk4(F,x0,t1,t2,n,p)

Args:
    F (function): Function that returns the derivatives of the state vector.
    x0 (array): Initial state vector.
    t1 (float): Initial time.
    t2 (float): Final time.
    n (int): Number of steps.
    p (array): Parameters.
"""
function rk4(F,x0,t1,t2,n,p)
    t = t1
    h = (t2-t1)/n
    x = x0

    for i in 1:n
        k1 = F(x,p,t)
        k2 = F(x .+ k1*(h/2),p,t+(h/2))
        k3 = F(x .+ k2*(h/2),p,t+(h/2))
        k4 = F(x .+ h*k3,p,t+h)

        x += h*(k1 .+ 2*k2 .+ 2*k3 .+ k4)/6
        t += h
    end
    return x
end

@doc raw"""Solves (system of) Ordinary differential equations using the Euler method.
    function euler(F,x0,t1,t2,n,p)

Args:
    F (function): Function that returns the derivatives of the state vector.
    x0 (array): Initial state vector.
    t1 (float): Initial time.
    t2 (float): Final time.
    n (int): Number of steps.
    p (array): Parameters."""
function euler(F,x0,t1,t2,n,p)
    t = t1
    h = (t2-t1)/n
    x = x0

    for i in 1:n
        x += h*F(x,p,t)
        t += h
    end
    return x
end

export rk4, euler