"""
    rotation_between(u, v)

Compute the quaternion that rotates vector `u` so that it aligns with vector
`v`, along the geodesic (shortest path).
"""
function rotation_between end

function rotation_between(u::StaticVector{2}, v::StaticVector{2})
    c = complex(v[1], v[2]) / complex(u[1], u[2])
    iszero(c) && throw(ArgumentError("Input vectors must be nonzero and finite."))
    isfinite(c) || throw(ArgumentError("Input vectors must be nonzero and finite."))
    theta = Base.angle(c)
    return Angle2d(theta)
end

function rotation_between(u::StaticVector{3}, v::StaticVector{3})
    # Robustified version of implementation from https://www.gamedev.net/topic/429507-finding-the-quaternion-betwee-two-vectors/#entry3856228
    normprod = sqrt(dot(u, u) * dot(v, v))
    T = typeof(normprod)
    normprod < eps(T) && throw(ArgumentError("Input vectors must be nonzero."))
    w = normprod + dot(u, v)
    v = abs(w) < 100 * eps(T) ? perpendicular_vector(u) : cross(u, v)
    @inbounds return QuatRotation(w, v[1], v[2], v[3]) # relies on normalization in constructor
end

function rotation_between(u::StaticVector{N}, v::StaticVector{N}) where N
    e1 = normalize(u)
    e2 = normalize(v - e1 * dot(e1, v))
    any(isnan, e1) && throw(ArgumentError("Input vectors must be nonzero and finite."))
    any(isnan, e2) && throw(ArgumentError("Input vectors must be nonzero and finite."))
    c = dot(e1, normalize(v))
    s = sqrt(1-c^2)
    P = [e1 e2 svd([e1 e2]'; full=true).Vt[StaticArrays.SUnitRange(3,N),:]']
    Q = one(MMatrix{N,N})
    Q[1,1] = c
    Q[1,2] = -s
    Q[2,1] = s
    Q[2,2] = c
    R = RotMatrix(P*Q*P')
    return R
end
