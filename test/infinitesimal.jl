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

    @testset "error case" begin
        for T in types_2d
            @test_throws BoundsError zero(InfinitesimalAngle2d)[5]
            @test_throws BoundsError zero(InfinitesimalAngle2d)[2,3]
            @test_throws BoundsError zero(InfinitesimalAngle2d)[3,1]
        end

        for T in types_3d
            @test_throws BoundsError zero(InfinitesimalAngle2d)[10]
            @test_throws BoundsError zero(InfinitesimalAngle2d)[2,4]
            @test_throws BoundsError zero(InfinitesimalAngle2d)[4,1]
        end

        @test_throws DimensionMismatch InfinitesimalAngle2d(1) + InfinitesimalRotationVec(2,3,4)
    end

    @testset "type promotion" begin
        for (T, N) in ((InfinitesimalAngle2d, 2), (InfinitesimalRotationVec, 3))
            R = InfinitesimalRotMatrix

            @test zero(T)      + zero(R{N})        isa R{N, Float64}
            @test zero(T{Int}) + zero(R{N,Int})    isa R{N, Int}
            @test zero(T{Int}) + zero(T{BigFloat}) isa T{BigFloat}

            @test zero(T)      - zero(R{N})        isa R{N, Float64}
            @test zero(T{Int}) - zero(R{N,Int})    isa R{N, Int}
            @test zero(T{Int}) - zero(T{BigFloat}) isa T{BigFloat}

            @test 42  * zero(T)      isa T{Float64}
            @test 42  * zero(T{Int}) isa T{Int}
            @test 1.2 * zero(T{Int}) isa T{Float64}
        end
    end
end
