# Rotation around an axis

In `Rotations.jl`, there's two types to represent a rotation with given an axis and an angle.

* `AngleAxis`
* `RotationVec`

## `AngleAxis`
A 3D rotation with fields `theta`, `axis_x`, `axis_y`, and `axis_z` to store the rotation angle and axis of the rotation.
Like all other types in this package, once it is constructed it acts and behaves as a 3×3 `AbstractMatrix`.
The axis will be automatically renormalized by the constructor to be a unit vector, so that `theta` always represents the rotation angle in radians.

## `RotationVec`
A 3D rotation encoded by an angle-axis representation as `angle * axis`.
This type is used in packages such as [OpenCV](http://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html#void%20Rodrigues%28InputArray%20src,%20OutputArray%20dst,%20OutputArray%20jacobian%29).

Note: If you're differentiating a Rodrigues Vector check the result is what you expect at theta = 0.
The first derivative of the rotation *should* behave, but higher-order derivatives of it (as well as parameterization conversions) should be tested.
The Stereographic Quaternion Projection (`MRP`) is the recommended three parameter format for differentiation.
