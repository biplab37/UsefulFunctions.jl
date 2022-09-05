## This file contains the code to write a bivariate function as a sum of separable functions.
## f(x,y) ≈ ∑ᵢ gᵢ(x)gᵢ(y)

function next_func(func::Function,p0)
    g(x,y) = func(x,y) - func(x,p0)*func(p0,y)/func(p0,p0)
    return g
end

function terms3(f::Function,p)

    F1(x,y) = f(x,y)
    F2 = next_func(F1,p[1])
    F3 = next_func(F2,p[2])

    f1(q) = F1(q,p[1])
    f2(q) = F2(q,p[2])
    f3(q) = F3(q,p[3])
    return f1,f2,f3
end

export terms3