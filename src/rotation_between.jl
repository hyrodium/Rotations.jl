"""
    rotation_between(from, to)

Compute the quaternion that rotates vector `from` so that it aligns with vector
`to`, along the geodesic (shortest path).
"""
function rotation_between(from::SVector{3}, to::SVector{3})
    # Robustified version of implementation from https://www.gamedev.net/topic/429507-finding-the-quaternion-betwee-two-vectors/#entry3856228
    normprod = sqrt(dot(from, from) * dot(to, to))
    T = typeof(normprod)
    normprod < eps(T) && throw(ArgumentError("Input vectors must be nonzero."))
    w = normprod + dot(from, to)
    v = abs(w) < 100 * eps(T) ? perpendicular_vector(from) : cross(from, to)
    @inbounds return QuatRotation(w, v[1], v[2], v[3]) # relies on normalization in constructor
end
