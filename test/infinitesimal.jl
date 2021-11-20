@testset "infinitesimal" begin
    all_types = (InfinitesimalRotMatrix{3}, InfinitesimalRotationVec,
                 InfinitesimalRotMatrix{2}, InfinitesimalAngle2d)
    types_2d = (InfinitesimalRotMatrix{2}, InfinitesimalAngle2d)
    types_3d = (InfinitesimalRotMatrix{3}, InfinitesimalRotationVec)

    @testset "zero" begin
        for T in all_types
            r = zero(T)
            @test r isa InfinitesimalRotation
            @test r === zero(r)

        end
    end

    @testset "one" begin
        for T in all_types
            r = one(T)
            @test r isa SMatrix
            @test r == one(r)
        end
    end

    @testset "minus" begin
        for T in all_types
            r = T(log(rand(typeof(exp(zero(T))))))
            @test r isa T
            @test -r isa T
            @test r' isa T
            @test transpose(r) isa T

            @test -r == -SMatrix(r)
            @test r' == SMatrix(r)'
            @test transpose(r) == transpose(SMatrix(r))

            @test r-r == r+(-r) == zero(r)
            @test -r == r' == transpose(r)
        end
    end

    @testset "multiply" begin
        for T in all_types
            r = T(log(rand(typeof(exp(zero(T))))))
            a = 4.2
            @test r isa T
            @test r*1 isa T
            @test 1*r isa T
            @test a*r isa T

            @test r*1 == 1*r == r
            @test r*2 == 2*r == r+r
            @test r*a == a*r == a*SMatrix(r)
        end
    end
end
