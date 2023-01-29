var documenterSearchIndex = {"docs":
[{"location":"diffeqn/#Algorithm-to-solve-Differential-equations","page":"Algorithm to solve Differential equations","title":"Algorithm to solve Differential equations","text":"","category":"section"},{"location":"diffeqn/#Euler-Method","page":"Algorithm to solve Differential equations","title":"Euler Method","text":"","category":"section"},{"location":"diffeqn/","page":"Algorithm to solve Differential equations","title":"Algorithm to solve Differential equations","text":"UsefulFunctions.euler","category":"page"},{"location":"diffeqn/#UsefulFunctions.euler","page":"Algorithm to solve Differential equations","title":"UsefulFunctions.euler","text":"Solves (system of) Ordinary differential equations using the Euler method.     function euler(F,x0,t1,t2,n,p)\n\nArgs:     F (function): Function that returns the derivatives of the state vector.     x0 (array): Initial state vector.     t1 (float): Initial time.     t2 (float): Final time.     n (int): Number of steps.     p (array): Parameters.\n\n\n\n\n\n","category":"function"},{"location":"diffeqn/#th-order-Runge-Kutta-Method","page":"Algorithm to solve Differential equations","title":"4th order Runge-Kutta Method","text":"","category":"section"},{"location":"diffeqn/","page":"Algorithm to solve Differential equations","title":"Algorithm to solve Differential equations","text":"UsefulFunctions.rk4","category":"page"},{"location":"diffeqn/#UsefulFunctions.rk4","page":"Algorithm to solve Differential equations","title":"UsefulFunctions.rk4","text":"Solves (system of) Ordinary differential equations using the Runge-Kutta method.     function rk4(F,x0,t1,t2,n,p)\n\nArgs:     F (function): Function that returns the derivatives of the state vector.     x0 (array): Initial state vector.     t1 (float): Initial time.     t2 (float): Final time.     n (int): Number of steps.     p (array): Parameters.\n\n\n\n\n\n","category":"function"},{"location":"indices/#List-of-Functions","page":"Indices","title":"List of Functions","text":"","category":"section"},{"location":"indices/","page":"Indices","title":"Indices","text":"List of the functions available","category":"page"},{"location":"indices/","page":"Indices","title":"Indices","text":"Order = [:function]","category":"page"},{"location":"root_finding/#Root-Finding-Algorithms","page":"Root Finding","title":"Root Finding Algorithms","text":"","category":"section"},{"location":"root_finding/#Bisection-Method","page":"Root Finding","title":"Bisection Method","text":"","category":"section"},{"location":"root_finding/","page":"Root Finding","title":"Root Finding","text":"UsefulFunctions.bisection","category":"page"},{"location":"root_finding/#UsefulFunctions.bisection","page":"Root Finding","title":"UsefulFunctions.bisection","text":"bisection(func::Function,start::Number,finish::Number,iteration::Integer=20)\n\nFinds the zero of a function func inside a given interval (start,finish). Might fail  if there are multiple zeros or no zeros inside the interval.\n\n\n\n\n\n","category":"function"},{"location":"root_finding/#Newton-Raphson","page":"Root Finding","title":"Newton-Raphson","text":"","category":"section"},{"location":"root_finding/","page":"Root Finding","title":"Root Finding","text":"NewtonRaphson","category":"page"},{"location":"root_finding/#UsefulFunctions.NewtonRaphson","page":"Root Finding","title":"UsefulFunctions.NewtonRaphson","text":"NewtonRaphson(f::Function,initial_guess;tol=1e-3,ϵ=1e-3,maxiter=20)\n\nUses Newton Raphson Method to find zero of a given function f starting from an initial_guess.\n\nx_n+1 = x_n - fracf(x_n)f(x_n)\n\nwhere discreet derivatives are used when functional form of the derivative are not given as input.\n\n\n\n\n\n","category":"function"},{"location":"#Useful-Functions","page":"Home","title":"Useful Functions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"A package containing some of the functions that I regularly use in my code.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Root Finding Algorithms\nSolving Differential Equations","category":"page"}]
}
