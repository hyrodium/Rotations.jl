@testset "Eigen" begin
    all_types = (RotMatrix{3}, AngleAxis, RotationVec,
                 UnitQuaternion, RodriguesParam, MRP,
                 RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                 RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
                 RotX, RotY, RotZ,
                 RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)

    for T in all_types, F in (one, rand)
        R = F(T)
        λs = eigvals(R)
        vs = eigvecs(R)
        E = eigen(R)
        v1 = vs[:,1]
        v2 = vs[:,2]
        v3 = vs[:,3]
        @test R * vs ≈ transpose(λs) .* vs
        @test norm(v1) ≈ 1
        @test norm(v2) ≈ 1
        @test norm(v3) ≈ 1
        @test E.values == λs
        @test E.vectors == vs
    end
end
