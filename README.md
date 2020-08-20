# This is a quick hack to compute accurate adjoint gradients for 3D photonic inverse design using MEEP FDTD.
The key ingredients in my hack are:
1. directly populate the Yee-grid through epsilon_func and amp_func <br/>
2. directly access the dft fields with yee_grid=True

# Test

The following is a 3D test structure consisting of two dielectric layers. Each layer has a randomized surface geometry. The degrees of freedom are the thickness of the dielectric at each point on each layer.

![alt text](https://github.com/zlin-opt/meep-adjoint-3d/blob/master/eps3d.png?raw=true)

The adjoint code is tested by sweeping the thickness of one of the dielectric layers at some random pixel point, and computing and comparing both the adjoint and finite-difference gradients. Note that in this example, I used a point current at some point above the dielectric and then computed the integrated field intensity |E|^2 at a plane beneath the dielectric. 

![alt text](https://github.com/zlin-opt/meep-adjoint-3d/blob/master/adj_diff.png?raw=true)

The error in percentage between the adjoint and finite-difference calculations.

![alt text](https://github.com/zlin-opt/meep-adjoint-3d/blob/master/adj_diff-err.png?raw=true)
