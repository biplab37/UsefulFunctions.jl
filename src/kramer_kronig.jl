## This file contains code related to the kramer-kronig transformation

function kramers_kronig_symmetric_shifted(func::Function, ω::Number, cutoff, integrate::Function)
    integrand(nu) = 2 * nu * func(nu) * (PrincipalValue(nu^2 - ω^2) - PrincipalValue(nu^2)) / π
    return integrate(integrand, 0.0, cutoff)
end

function kramers_kronig_asymmetric_shifted(func::Function, ω::Number, cutoff_below, cutoff_above, integrate::Function)
    integrand(nu) = func(nu) * (PrincipalValue(nu - ω) - PrincipalValue(nu)) / π
    return integrate(integrand, cutoff_below, cutoff_above)
end
