module UsefulFunctions

"""
	DiracDelta(input::Float64,δ::Float64 = 1e-3)

Dirac Delta Function that can be used in a numerical integration.
"""

function DiracDelta(input, δ = 1e-3)
	return exp(-input^2/(4*δ))/(2*sqrt(pi*δ))
end

"""
	PrincipleValue(x::Float64,ϵ::Float64 = 1e-3)

Useful function in Principle Value Integration.
"""
function PrincipleValue(x,ϵ = 1e-3)
	if abs(x)<ϵ
		return 0.0
	else
		return 1/x
	end
end

"""
	fzero(func::Function, start::Float64, finish::Float64, iteration::Int64=20)

Returns the value of the argument between `start` and `finish` where the given function `func` becomes zero using midpoint method.
"""
function fzero(func, start, finish, iteration=20)
	mid = (start + finish)/2.0

    for i in 1:iteration
        if func(mid)*func(start) > 0
            start = mid
        else
        	finish = mid
        end
    	mid = (start + finish)/2.0
    end
    return mid
end

"""
	numberF(temp::Float64, μ::Float64, Energy::Float64)

Returns the number density at a given `Energy`, chemical potential `μ` and `temp` for Fermionic Species.
"""
function numberF(temp, μ, Energy)
	β = 1/temp
    return 1/(1+exp(β*(Energy - μ)))
end

"""
	numberB(temp::Float64, μ::Float64, Energy::Float64)

Returns the number density at a given `Energy`, chemical potential `μ` and `temp` for Bosonic Species.
"""
function numberB(temp, μ, Energy)
	β = 1/temp
    return 1/(1 - exp(β*(Energy - μ)))
end

function itp(momentum::Float64,n::Int64,list::Array)

    index1 = Int64(floor(n*momentum))
    index2 = Int64(ceil(n*momentum))

    if index1==0
        return list[1]
    elseif index2>n
        return list[n]
    else
        if index1 == index2
            return list[index1]
        else
            return list[index1] + (list[index2] - list[index1])*(n*momentum - index1)
        end
    end
end

@doc raw"""
    function interp(n::Int64,list::Array)
This function returns the interpolated function between 0 and 1 given an Array.
"""
function interp(list::Array)
	n = length(list)
    itp2(momentum::Float64) = itp(momentum,n,list)
    return itp2
end

include("root_finding.jl")
include("diffeqn.jl")
include("separable.jl")

export DiracDelta, PrincipleValue, numberF, fzero, interp

end # module
