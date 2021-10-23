@testset "Eigen" begin
    all_types = (RotMatrix{3}, AngleAxis, RotationVec,
                 UnitQuaternion, RodriguesParam, MRP,
                 RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                 RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ,
                 RotX, RotY, RotZ,
                 RotXY, RotYZ, RotZX, RotXZ, RotYX, RotZY)

    for type in all_types
        R = rand(type)
        λs = eigvals(R)
        vs = eigvecs(R)
        @test R * vs ≈ transpose(λs) .* vs
    end
end
