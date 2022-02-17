#=
The following tests check the distribution of `rand(RotType)`.
To visualize the distribution, try the following script.

    using Rotations, Plots
    RotType = RotMatrix3{Float64}

    # Sampling
    n = 10000
    rs = rand(RotType,n)
    q = rand(RotType)
    angles = rotation_angle.([r/q for r in rs])/π

    # Cumulative distribution function
    cdf_sampled(t) = count(<(t), angles)/n
    cdf_true(t) = t-sinpi(t)/π

    # Probability density function and histogram of sampling distribution
    histogram(angles, normed=true, label=false)
    plot!(t->1-cospi(t), 0, 1, label=false)

    # Cumulative distribution function
    plot(cdf_sampled, 0, 1, label=false)
    plot!(cdf_true, 0, 1, label=false)
=#

@testset "Distribution" begin

    # TODO: consider to remove Euler rotations (#155)
    Type_SO3 =  (RotMatrix3, RotMatrix{3}, AngleAxis, RotationVec,
                QuatRotation, MRP, RodriguesParam,
                RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ)

    Type_SO2 = (RotMatrix2, RotMatrix{2}, Angle2d, RotX, RotY, RotZ)

    # Number of sampling
    N = 100000

    @testset "Distribution $(RotType)" for RotType in Type_SO3
        # Sampling
        rs = rand(RotType, N)
        q = rand(RotType)
        angles = rotation_angle.([r/q for r in rs])/π

        # Check sampled rotations are rotations
        @test all(isrotation.(rs))
        @test isrotation(q)

        # Cumulative distribution function
        cdf_sampled(t) = count(<(t), angles)/N
        cdf_true(t) = t-sinpi(t)/π

        # Check the CDFs are approximately equal
        ts = 0:0.01:1
        @test norm(cdf_true.(ts) - cdf_sampled.(ts), Inf) < 0.01
    end

    @testset "distribution $(RotType)" for RotType in Type_SO2
        # Sampling
        rs = rand(RotType, N)
        q = rand(RotType)
        angles = rotation_angle.([r/q for r in rs])/2π
        angles = mod.(angles, 1)

        # Check sampled rotations are rotations
        @test all(isrotation.(rs))
        @test isrotation(q)

        # Cumulative distribution function
        cdf_sampled(t) = count(<(t), angles)/N
        cdf_true(t) = t

        # Check the CDFs are approximately equal
        ts = 0:0.01:1
        @test norm(cdf_true.(ts) - cdf_sampled.(ts), Inf) < 0.01
    end

    @testset "distribution $(N)-dimensinal RotMatrix" for N in 2:5
        # Sampling
        RotType = RotMatrix{N}
        n = 1000000
        rs = rand(RotType, n)
        q = rand(RotType)
        norms1 = norm.([(r/q - one(SMatrix{N,N})) for r in rs])
        norms2 = norm.([(r - one(SMatrix{N,N})) for r in rs])

        # Check sampled rotations are rotations
        @test all(isrotation.(rs))
        @test isrotation(q)

        # Cumulative distribution function
        cdf_sampled1(t) = count(<(t), norms1)/n
        cdf_sampled2(t) = count(<(t), norms2)/n

        # Check the CDFs are approximately equal
        ts = 0:0.01:1
        @test norm(cdf_sampled1.(ts) - cdf_sampled2.(ts), Inf) < 0.01
    end
end
