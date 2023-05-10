"""
    rotation_between(u, v)

Compute the quaternion that rotates vector `u` so that it aligns with vector
`v`, along the geodesic (shortest path).
"""
function rotation_between(u::SVector{3}, v::SVector{3})
    # Robustified version of implementation from https://www.gamedev.net/topic/429507-finding-the-quaternion-betwee-two-vectors/#entry3856228
    normprod = sqrt(dot(u, u) * dot(v, v))
    T = typeof(normprod)
    normprod < eps(T) && throw(ArgumentError("Input vectors must be nonzero."))
    w = normprod + dot(u, v)
    v = abs(w) < 100 * eps(T) ? perpendicular_vector(u) : cross(u, v)
    @inbounds return QuatRotation(w, v[1], v[2], v[3]) # relies on normalization in constructor
end
