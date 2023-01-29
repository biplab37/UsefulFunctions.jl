using UsefulFunctions, Test

@testset "Integration" begin
    @testset "Trapizoid" begin
        @test trapizoid(x->x^2, 0., 1., 10000) ≈ 0.3333333333333333
    end
    @testset "Simpson" begin
        @test simpson(x->x^2, 0., 1., 1000) ≈ 0.3333333333333333
        @test simpson(x->x^2, 0., 1., 10000) ≈ 0.3333333333333333
    end
end