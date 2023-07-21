# gp_only_sphere_test
Approximate volume of the unit sphere using gp-only integration

Approximates the volume of the unit sphere using gradient path only integration.
Uses six gradient paths which are the positive and negative x, y, and z axes.

getpaths(n,r) generates the six paths with n points starting on an interior sphere of radius r.

int_vol(n,r) approximates volume of the unit sphere using n points on the gradient paths for all integration. Gradient paths begin on an interior sphere of radius r. Returns relative error.

errortable() generates a table of error and convergence results varying the choice of n. Uses r=.05.

