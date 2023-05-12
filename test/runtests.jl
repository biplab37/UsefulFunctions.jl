using UsefulFunctions, Test

sq(x) = x^2

f(x) = x^2 - 1
df(x) = 2x



x_init = 0.0
x_final = 2.0

x0 = 0.1

@testset "Integration" begin
    @testset "Trapizoid" begin
        @test trapizoid(sq, 0.0, 1.0, 100) ≈ 1 / 3 atol = 1e-4
        @test trapizoid(sq, 0.0, 1.0, 1000) ≈ 1 / 3 atol = 1e-5
        @test trapizoid(sq, 0.0, 1.0, 10000) ≈ 1 / 3
    end
    @testset "Simpson" begin
        @test simpson(sq, 0.0, 1.0, 100) ≈ 0.3333333333333333 atol = 1e-5
        @test simpson(sq, 0.0, 1.0, 1000) ≈ 0.3333333333333333
        @test simpson(sq, 0.0, 1.0, 10000) ≈ 0.3333333333333333
    end
end

@testset "Root Finding" begin
    @testset "Bisection" begin
        @test bisection(f, 0.0, 2.0) ≈ 1.0 atol = 1e-5
        @test bisection(sq, 0.0, 2.0) ≈ 0.0 atol = 1e-5
        @test bisection(sq, 0.0, 2.0, 15) == 2^(-15)
    end
    @testset "fzero" begin
        @test fzero(f, 0.0, 2.0) ≈ 1.0 atol = 1e-5
        @test fzero(sq, 0.0, 2.0) ≈ 0.0 atol = 1e-5
        @test fzero(sq, 0.0, 2.0, 15) == 2^(-15)
    end
    @testset "NewtonRaphson" begin
        @test NewtonRaphson(f, x0) ≈ 1.0 atol = 1e-5
    end
    @testset "Newton Method" begin
        @test NewtonMethod(f, df, x0) ≈ 1.0
        @test NewtonMethod(f, df, x0, 10) ≈ 1.0
        @test typeof(NewtonMethod(f, df, x0)) == Float64
        @test typeof(NewtonMethod(f, df, x0 * (1 + im))) == Complex{Float64}
        @test abs(real(NewtonMethod(f, df, x0 * (1 + im)))) ≈ 1.0 atol = 1e-6
        @test abs(imag(NewtonMethod(f, df, x0 * (1 + im)))) ≈ 0.0
        @test NewtonMethod(x -> [x[1]^2 - 2, x[2]^2 - 4], x -> [2x[1] 0; 0 2x[2]], [x0, x0]) ≈ [√2, 2]
    end
    @testset "Brent" begin
        @test brent(f, 0.0, 2.0) ≈ 1.0 atol = 1e-5
        @test brent(sq, 0.0, 2.0) ≈ 0.0 atol = 1e-5
    end
    @testset "Broyden" begin
        @test broyden(f, df, x0) ≈ 1.0 atol = 1e-5
        @test broyden(x -> [x[1]^2 - 2, x[2]^2 - 4], x -> [2x[1] 0; 0 2x[2]], [1.0, 3.0]) ≈ [√2, 2]
    end
end

@testset "Number Density" begin
    @test numberF(0.001, 0.0, 0.1) ≈ 0.0 atol = 1e-10
    @test numberF(0.001, 0.2, 0.1) ≈ 1.0 atol = 1e-10
    @test numberB(0.001, 0.0, 0.1) ≈ 0.0 atol = 1e-10
    @test numberB(0.001, 0.2, 0.1) ≈ -1.0 atol = 1e-10
end

F(x, p, t) = p * x
t1 = 0.0
t2 = 1.0
euler(F, x0, t1, t2, 500, 0.01)

@testset "Differential Eqn" begin
    @testset "Euler" begin
        @test euler(F, x0, t1, t2, 50, 0.01) ≈ 0.10 atol = 1e-2
    end
    @testset "RK4" begin
        @test rk4(F, x0, t1, t2, 50, 0.01) ≈ 0.10 atol = 1e-2
    end
end

@testset "General Functions" begin
    @testset "Dirac Delta" begin
        @test DiracDelta(0.5) ≈ 0.0 atol = 1e-10
        @test DiracDelta(0.0) > 1.0
        @test DiracDelta(0.1, 1e-4) > 1e-10
    end

    @testset "PrincipalValue" begin
        @test PrincipalValue(0.0) == 0.0
        @test PrincipalValue(0.01, 0.005) == 100.0
        @test PrincipalValue(0.001, 0.01) == 0.0
    end

end

@testset "Interpolation" begin
    @test interp(collect(1:10))(0.05) == 1.0
    @test interp(collect(1:10))(0.95) == 9.5
    @test interp(collect(1:10))(0.45) == 4.5
end
