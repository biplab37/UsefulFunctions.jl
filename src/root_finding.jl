function NewtonMethod(func::Function, Jacobian,initial_guess,iteration=1000)
	for i in 1:iteration
		initial_guess -= inv(Jacobian(initial_guess))*func(initial_guess)
	end

	return initial_guess
end

function NewtonRaphson(f::Function,initial_guess;tol=1e-3,ϵ=1e-3,maxiter=20)
	x_old = initial_guess
	x_new = x_old
	error = 1
	for _ in 1:maxiter
		fx = f(x_old)
		x_new = x_old - fx*ϵ/(f(x_old + ϵ) - fx)
		if abs(x_new - x_old) < tol
			return x_new
		end
		x_old = x_new
		x_new = nothing
	end

	return x_new
end

function brent(f::Function, x0::Number, x1::Number, args::Tuple=(); xtol::AbstractFloat=1e-4, ytol=1e-6, maxiter::Integer=20)
	EPS = eps(Float64)
	y0 = f(x0,args...)
	y1 = f(x1,args...)
	if abs(y0) < abs(y1)
		# Swap lower and upper bounds.
		x0, x1 = x1, x0
		y0, y1 = y1, y0
	end
	x2 = x0
	y2 = y0
	x3 = x2
	bisection = true
	for _ in 1:maxiter
		# x-tolerance.
		if abs(x1-x0) < xtol
			return x1
		end

		# Use inverse quadratic interpolation if f(x0)!=f(x1)!=f(x2)
		# and linear interpolation (secant method) otherwise.
		if abs(y0-y2) > ytol && abs(y1-y2) > ytol
			x = x0*y1*y2/((y0-y1)*(y0-y2)) + x1*y0*y2/((y1-y0)*(y1-y2)) + x2*y0*y1/((y2-y0)*(y2-y1))
		else
			x = x1 - y1 * (x1-x0)/(y1-y0)
		end

		# Use bisection method if satisfies the conditions.
		delta = abs(2EPS*abs(x1))
		min1 = abs(x-x1)
		min2 = abs(x1-x2)
		min3 = abs(x2-x3)
		if (x < (3x0+x1)/4 && x > x1) ||	(bisection && min1 >= min2/2) ||	(!bisection && min1 >= min3/2) ||	(bisection && min2 < delta) ||	(!bisection && min3 < delta)
			x = (x0+x1)/2
			bisection = true
		else
			bisection = false
		end

		y = f(x,args...)
		# y-tolerance.
		if abs(y) < ytol
			return x
		end
		x3 = x2
		x2 = x1
		if sign(y0) != sign(y)
			x1 = x
			y1 = y
		else
			x0 = x
			y0 = y
		end
		if abs(y0) < abs(y1)
			# Swap lower and upper bounds.
			x0, x1 = x1, x0
			y0, y1 = y1, y0
		end
	end
	error("Max iteration exceeded")
end

function broyden(fun, jaco, x, iter=50, ftol=1e-7, verbose=false)
    msg = "Maximum number of iterations reached."
    J = jaco(x)
    for cont = 1:iter
        if det(J) < ftol
            x = nothing
            msg = "Derivative near to zero."
            break
        end
        if verbose
            println("n: $(cont), x: $(x)")
        end
        f_old = fun(x)
        dx = -J\f_old
        x = x + dx
        f = fun(x)
        df = f - f_old
        J = J + (df - J*dx) * dx'/ (dx' * dx)
        if norm(f) < ftol
            msg = "Root found with desired accuracy."
            break
        end
    end
    return x
end

function bisection(func, start, finish, iteration=20)
	mid = (start + finish)/2.0

    for i in 1:iteration
        if func(mid)*func(finish) > 0
            finish = mid
        else
        	start = mid
        end
    	mid = (start + finish)/2.0
    end
    return mid
end

export bisection, NewtonMethod, NewtonRaphson, brent, broyden