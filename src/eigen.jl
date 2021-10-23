function LinearAlgebra.eigvals(R::Rotation{3})
    θ = rotation_angle(R)
    return SVector(exp(θ*im), exp(-θ*im), 1)
end

function LinearAlgebra.eigvecs(R::Rotation{3})
    n3 = normalize(rotation_axis(R))
    n1 = normalize(perpendicular_vector(n3))
    n2 = normalize(n3×n1)
    v1 = n1*im + n2
    v2 = n1 + n2*im
    v3 = n3
    return hcat(v1,v2,v3)
end
