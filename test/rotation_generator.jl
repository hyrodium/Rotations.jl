@testset "rotation generator" begin
    all_types = (RotMatrixGenerator{3}, RotationVecGenerator,
                 RotMatrixGenerator{2}, Angle2dGenerator)
    types_2d = (RotMatrixGenerator{2}, Angle2dGenerator)
    types_3d = (RotMatrixGenerator{3}, RotationVecGenerator)

    @testset "constructor-2d" begin
        m = rand(2,2)
        s1 = RotMatrixGenerator{2}(m - m')
        s2 = RotMatrixGenerator{2, BigFloat}(m - m')
        s3 = RotMatrixGenerator{2}(0)
        s4 = RotMatrixGenerator{2, BigFloat}(0)
        s5 = RotMatrixGenerator(0)
        s6 = RotMatrixGenerator(BigFloat(0))
        s7 = Angle2dGenerator(0)
        s8 = Angle2dGenerator(BigFloat(0))
        s9 = Angle2dGenerator{BigFloat}(0)
        @test s1 isa RotMatrixGenerator{2, Float64}
        @test s2 isa RotMatrixGenerator{2, BigFloat}
        @test s3 isa RotMatrixGenerator{2, Int}
        @test s4 isa RotMatrixGenerator{2, BigFloat}
        @test s5 isa RotMatrixGenerator{2, Int}
        @test s6 isa RotMatrixGenerator{2, BigFloat}
        @test s7 isa Angle2dGenerator{Int}
        @test s8 isa Angle2dGenerator{BigFloat}
        @test s9 isa Angle2dGenerator{BigFloat}
    end

    @testset "constructor-3d" begin
        m = rand(3,3)
        s1 = RotMatrixGenerator{3}(m - m')
        s2 = RotMatrixGenerator{3, BigFloat}(m - m')
        @test s1 isa RotMatrixGenerator{3, Float64}
        @test s2 isa RotMatrixGenerator{3, BigFloat}
    end

    @testset "zero" begin
        for T in all_types
            r = zero(T)
            @test r isa RotationGenerator
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

        @test one(RotationGenerator{2}) isa SMatrix{2, 2, Float64}
        @test one(RotationGenerator{2,BigFloat}) isa SMatrix{2, 2, BigFloat}
        @test one(RotationGenerator{3}) isa SMatrix{3, 3, Float64}
        @test one(RotationGenerator{3,BigFloat}) isa SMatrix{3, 3, BigFloat}
    end

    @testset "minus" begin
        for T in all_types
            # TODO: These should be replaced with `r = rand(T)`
            if T in types_2d
                r = T(Angle2dGenerator(1.2))
            elseif T in types_3d
                r = T(RotationVecGenerator(1.2,-0.8,0.1))
            end
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
            # TODO: These should be replaced with `r = rand(T)`
            if T in types_2d
                r = T(Angle2dGenerator(1.2))
            elseif T in types_3d
                r = T(RotationVecGenerator(1.2,-0.8,0.1))
            end
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
            # TODO: These should be replaced with `r = rand(T)`
            if T in types_2d
                r = T(Angle2dGenerator(1.2))
            elseif T in types_3d
                r = T(RotationVecGenerator(1.2,-0.8,0.1))
            end
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
            # TODO: These should be replaced with `r = rand(T)`
            if T in types_2d
                r = T(Angle2dGenerator(1.2))
            elseif T in types_3d
                r = T(RotationVecGenerator(1.2,-0.8,0.1))
            end
            @test r isa T
            @test r*r isa SMatrix
            @test r/r isa SMatrix
            @test r^2 isa SMatrix
        end
    end

    @testset "error case" begin
        for T in types_2d
            @test_throws BoundsError zero(Angle2dGenerator)[5]
            @test_throws BoundsError zero(Angle2dGenerator)[2,3]
            @test_throws BoundsError zero(Angle2dGenerator)[3,1]
        end

        for T in types_3d
            @test_throws BoundsError zero(Angle2dGenerator)[10]
            @test_throws BoundsError zero(Angle2dGenerator)[2,4]
            @test_throws BoundsError zero(Angle2dGenerator)[4,1]
        end

        @test_throws ErrorException zero(RotationGenerator)
        @test_throws ErrorException one(RotationGenerator)
        @test_throws ErrorException zero(RotMatrixGenerator)
        @test_throws ErrorException one(RotMatrixGenerator)

        @test_throws DimensionMismatch Angle2dGenerator(1) + RotationVecGenerator(2,3,4)
    end

    @testset "params" begin
        @test Rotations.params(Angle2dGenerator(1)) == [1]
        @test Rotations.params(RotationVecGenerator(2,3,4)) == [2,3,4]
    end

    @testset "type promotion" begin
        for (T, N) in ((Angle2dGenerator, 2), (RotationVecGenerator, 3))
            R = RotMatrixGenerator

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
        r = zero(RotMatrixGenerator{2})
        show(io, MIME("text/plain"), r)
        str = String(take!(io))
        @test startswith(str, "2×2 RotMatrixGenerator2{Float64}")

        rvec = RotationVecGenerator(1.0, 2.0, 3.0)
        show(io, MIME("text/plain"), rvec)
        str = String(take!(io))
        @test startswith(str, "3×3 RotationVecGenerator{Float64}") && occursin("(1.0, 2.0, 3.0):", str)
    end
end
