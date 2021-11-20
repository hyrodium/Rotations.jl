@testset "infinitesimal" begin
    all_types = (InfinitesimalRotMatrix{3}, InfinitesimalRotationVec,
                 InfinitesimalRotMatrix{2}, InfinitesimalAngle2d)
    types_2d = (InfinitesimalRotMatrix{2}, InfinitesimalAngle2d)
    types_3d = (InfinitesimalRotMatrix{3}, InfinitesimalRotationVec)

    @testset "constructor" begin
        m = rand(3,3)
        s1 = InfinitesimalRotMatrix{3}(m - m')
        s2 = InfinitesimalRotMatrix{3, BigFloat}(m - m')
        @test s1 isa InfinitesimalRotMatrix{3, Float64}
        @test s2 isa InfinitesimalRotMatrix{3, BigFloat}
    end

    @testset "zero" begin
        for T in all_types
            r = zero(T)
            @test r isa InfinitesimalRotation
            @test r === zero(r)
            @test zeros(T,2,3) == [zero(T) for i in 1:2, j in 1:3]
            @test zeros(T{BigFloat},2,3) == [zero(T{BigFloat}) for i in 1:2, j in 1:3]
            @test zeros(T,(2,3)) == [zero(T) for i in 1:2, j in 1:3]
            @test zeros(T{BigFloat},(2,3)) == [zero(T{BigFloat}) for i in 1:2, j in 1:3]
        end
    end

    @testset "one" begin
        for T in all_types
            r = one(T)
            @test r isa SMatrix
            @test r == one(r)
            @test ones(T,2,3) == [one(T) for i in 1:2, j in 1:3]
            @test ones(T{BigFloat},2,3) == [one(T{BigFloat}) for i in 1:2, j in 1:3]
            @test ones(T,(2,3)) == [one(T) for i in 1:2, j in 1:3]
            @test ones(T{BigFloat},(2,3)) == [one(T{BigFloat}) for i in 1:2, j in 1:3]
        end

        for T in types_2d
            r = one(T{BigFloat})
            @test r isa SMatrix{2,2,BigFloat}
        end

        for T in types_3d
            r = one(T{BigFloat})
            @test r isa SMatrix{3,3,BigFloat}
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

    @testset "division" begin
        for T in all_types
            r = T(log(rand(typeof(exp(zero(T))))))
            a = 4.2
            @test r isa T
            @test r/1 isa T
            @test r/a isa T

            @test r/1 == r
            @test r/2 + r/2 ≈ r
            @test r/a == SMatrix(r)/a
        end
    end

    @testset "matrix multiplication" begin
        for T in all_types
            r = T(log(rand(typeof(exp(zero(T))))))
            @test r isa T
            @test r*r isa SMatrix
            @test r/r isa SMatrix
            @test r^2 isa SMatrix
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

        @test_throws ErrorException zero(InfinitesimalRotation)
        @test_throws ErrorException one(InfinitesimalRotation)

        @test_throws DimensionMismatch InfinitesimalAngle2d(1) + InfinitesimalRotationVec(2,3,4)
    end

    @testset "type promotion" begin
        for (T, N) in ((InfinitesimalAngle2d, 2), (InfinitesimalRotationVec, 3))
            R = InfinitesimalRotMatrix

            @test zero(T)           + zero(R{N})        isa R{N, Float64}
            @test zero(R{N})        + zero(T)           isa R{N, Float64}
            @test zero(T{Int})      + zero(R{N,Int})    isa R{N, Int}
            @test zero(R{N,Int})    + zero(T{Int})      isa R{N, Int}
            @test zero(T{Int})      + zero(T{BigFloat}) isa T{BigFloat}
            @test zero(T{BigFloat}) + zero(T{Int})      isa T{BigFloat}

            @test zero(T)           - zero(R{N})        isa R{N, Float64}
            @test zero(R{N})        - zero(T)           isa R{N, Float64}
            @test zero(T{Int})      - zero(R{N,Int})    isa R{N, Int}
            @test zero(R{N,Int})    - zero(T{Int})      isa R{N, Int}
            @test zero(T{Int})      - zero(T{BigFloat}) isa T{BigFloat}
            @test zero(T{BigFloat}) - zero(T{Int})      isa T{BigFloat}

            @test 42  * zero(T)      isa T{Float64}
            @test 42  * zero(T{Int}) isa T{Int}
            @test 1.2 * zero(T{Int}) isa T{Float64}
        end
    end

    @testset "Testing show" begin
        io = IOBuffer()
        r = zero(InfinitesimalRotMatrix{2})
        show(io, MIME("text/plain"), r)
        str = String(take!(io))
        if VERSION ≥ v"1.6"
            @test startswith(str, "2×2 InfinitesimalRotMatrix2{Float64}")
        else
            @test startswith(str, "2×2 InfinitesimalRotMatrix{2,Float64,4}")
        end

        rvec = InfinitesimalRotationVec(1.0, 2.0, 3.0)
        show(io, MIME("text/plain"), rvec)
        str = String(take!(io))
        @test startswith(str, "3×3 InfinitesimalRotationVec{Float64}") && occursin("(1.0, 2.0, 3.0):", str)
    end
end
