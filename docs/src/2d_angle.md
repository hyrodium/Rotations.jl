# 2D Rotation by Angle

`Angle2d` is a type that internally has an angle as a parameter, and can be thought of as a lazy version of `RotMatrix{2}`.

```@setup angle2d
using Rotations
```

## `Angle2d`

**example**

```@repl angle2d
t = 1.2  # rotation angle
m = [cos(t) -sin(t);sin(t) cos(t)]
Angle2d(t)  # construct from matrix
RotMatrix{2}(t)
dump(Angle2d(t))
dump(RotMatrix{2}(t))
```
