# 3D Rotation around an axis

In `Rotations.jl`, there are two concrete types to represent a rotation with given an axis and an angle.

* `AngleAxis`
* `RotationVec`

## `AngleAxis`
A 3D rotation with fields `theta`, `axis_x`, `axis_y`, and `axis_z` to store the rotation angle and axis of the rotation.
Like all other types in this package, once it is constructed it acts and behaves as a 3×3 `AbstractMatrix`.
The axis will be automatically renormalized by the constructor to be a unit vector, so that `theta` always represents the rotation angle in radians.

```math
\begin{aligned}
R_{\bm{n}} (\theta)
&= \exp(\theta K(\bm{n})) \\
&= I + (\sin\theta) K(\bm{n})  +  (1-\cos\theta) K^2(\bm{n}) \\
&= \begin{pmatrix}
\cos\theta + n_x^2 (1-\cos\theta) & n_xn_y (1-\cos\theta) - n_z \sin\theta & n_zn_x (1-\cos\theta) + n_y \sin\theta \\
n_xn_y (1-\cos\theta) + n_z \sin\theta & \cos\theta + n_y^2 (1-\cos\theta) & n_yn_z (1-\cos\theta) - n_x \sin\theta \\
n_zn_x (1-\cos\theta) - n_y \sin\theta & n_yn_z (1-\cos\theta) + n_x \sin\theta & \cos\theta + n_z^2 (1-\cos\theta)
\end{pmatrix} \\
%
K(\bm{n})
&= \begin{pmatrix}
0    & -n_z & n_y \\
n_z  & 0    & -n_x \\
-n_y & n_x  & 0
\end{pmatrix} \\
%
{\bm{n}}
&= \begin{pmatrix}
n_x \\
n_y \\
n_z
\end{pmatrix}
%
\quad \left(\|\bm{n}\|= 1\right)
\end{aligned}
```

### Example

```@repl
using Rotations, LinearAlgebra
# 1/3 Rotation around (1/√3, 1/√3, 1/√3) vector
R = AngleAxis(2π/3, 1/√3, 1/√3, 1/√3)
# This matrix swaps the xyz coordinates
R * [1,2,3]
R^2
R^3
# These matrices are approximately equal
R^3 ≈ I(3)
```

## `RotationVec`
A 3D rotation encoded by an angle-axis representation as `angle * axis`.
This type is used in packages such as [OpenCV](http://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html#void%20Rodrigues%28InputArray%20src,%20OutputArray%20dst,%20OutputArray%20jacobian%29).

Note: If you're differentiating a Rodrigues Vector check the result is what you expect at theta = 0.
The first derivative of the rotation *should* behave, but higher-order derivatives of it (as well as parameterization conversions) should be tested.
The Stereographic Quaternion Projection (`MRP`) is the recommended three parameter format for differentiation.
