# Euler Angles

[Euler angles](https://en.wikipedia.org/wiki/Euler_angles) are a way to represent a rotation matrix with three rotations around cardinal axes.

In Rotations.jl, there are 12 concrete types for Euler angles.

* Proper Euler angles
    * `RotZXZ`, `RotXYX`, `RotYZY`, `RotZYZ`, `RotXZX`, `RotYXY`
* Tait–Bryan angles
    * `RotXYZ`, `RotYZX`, `RotZXY`, `RotXZY`, `RotZYX`, `RotYXZ`

In addition, Rotations.jl provides types that represent rotations in one or two axes.

* one axis
    * `RotX`, `RotY`, `RotZ`
* two axes
    * `RotXY`, `RotYZ`, `RotZX`, `RotXZ`, `RotZY`, `RotYX`

## Rotation around one axis

### Example
```@setup one_axis
using Rotations
```

Here's an example for `RotZ`:
```@repl one_axis
α = 1.2  # Rotation angle
R = RotZ(α)
Q = [cos(α) -sin(α) 0
      sin(α)  cos(α) 0
      0       0      1]
R == Q  # These matrices are equal
```

And more examples for `RotX` and `RotY`:

```@repl one_axis
RotX(α) == [1 0       0
            0 cos(α) -sin(α)
            0 sin(α)  cos(α)]
RotY(α) == [cos(α) 0 sin(α)
            0      1 0
           -sin(α) 0 cos(α)]
```

## Rotation around two axes

### Example
```@setup two_axis
using Rotations
```

```@repl two_axis
α, β = 1.2, 4.7  # Rotation angles
RotX(α) * RotY(β)
RotXY(α, β)
RotX(α) * RotY(β) == RotXY(α, β)
```

## Rotation around three axes

### Example
```@setup three_axis
using Rotations
```

```@repl three_axis
α, β, γ = 1.2, 4.7, -0.4  # Rotation angles
RotYXZ(α, β, γ)
RotY(α)*RotX(β)*RotZ(γ)
RotYXZ(α, β, γ) == RotY(α)*RotX(β)*RotZ(γ)
```
