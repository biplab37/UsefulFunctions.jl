module UsefulFunctions

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
	σ1(temp,μ, κ=0.05, M=1.0)

"""
function σ1(temp,μ, κ=0.05, M=1.0)
    β = 1/temp
	
    σ₁(σ) = (1 + 2*cosh(β*μ)*exp(-β*σ) + exp(-2*β*σ) - exp(β*(M - σ + (π*κ/σ))))
	
    result = fzero(σ₁,0.001,2)
	
    ## For the case when there is no zero in the interval
	if (result == 2.0 || result == 0.001)
        result = 0.0
    end
	
    return result
end

export PrincipleValue, numberF, σ1, fzero

end # module
