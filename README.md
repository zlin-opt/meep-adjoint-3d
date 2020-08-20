# This is a quick hack to compute accurate adjoint gradients for 3D photonic inverse design using MEEP FDTD.
The key ingredients in my hack are:
1. directly populate the Yee-grid through epsilon_func and amp_func <br/>
2. directly access the dft fields with yee_grid=True
