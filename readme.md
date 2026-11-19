# Algebraic curve reconstruction

Recover a binary shape from a coarse grid of measurements when its boundary is the zero set of a polynomial in two variables. The Bernstein method uses B-spline samples to estimate image moments in a triangular Bernstein basis, then fits the polynomial whose zero contour gives the boundary.

## Run the comparison

From the repository root, open MATLAB or GNU Octave. In Octave, load the `image`, `optim`, and `statistics` packages first. The corresponding MATLAB toolboxes are Image Processing Toolbox, Optimization Toolbox, and Statistics and Machine Learning Toolbox. The full comparison has been tested in Octave.

```matlab
% GNU Octave only
pkg load image
pkg load optim
pkg load statistics
```

```matlab
addpath(genpath(pwd))
rng(20261119)
compare_methods
mean(d_array, 1)
```

The experiment processes all 100 bundled 501 × 501 masks. For each mask it simulates coarse measurements, adds Gaussian noise at a target signal-to-noise ratio of 50 dB, and fits a degree-four boundary. The fixed seed makes the noise repeatable within the tested environment. The power, Bernstein, and Vetterli fits draw noise separately. Both Vetterli outputs use the same noisy measurements.

The boxplot shows how much each reconstruction overlaps the original mask. A Dice score of 1 means exact overlap and 0 means no overlap. The script also leaves a 100 × 4 matrix named `d_array` in the workspace, with one row per image. Its columns are power basis, triangular Bernstein basis, Vetterli quadratic-program fit, and the same Vetterli fit after refinement.

## Explore the code

- [Comparison script](laboratory/scripts/compare_methods.m) — experiment settings, input data, and scoring.
- [Triangular Bernstein reconstruction](bases/non_separable/functions/non_separable_reconstruction.m) — sampling, moment equations, and coefficient fitting.
- [Basis reproduction check](bases/non_separable/scripts/check_triangular_basis_reproduction.m) — a short numerical check of the Bernstein expansion using the spline kernel.
- [Bundled comparison images](data/statistics/images_02.mat) — the reference masks used by the comparison.

## Paper

Costanza Conti, Mariantonia Cotronei, Demetrio Labate, and Wilfredo Molina (2021). “Stable recovery of planar regions with algebraic boundaries in Bernstein form.” *Advances in Computational Mathematics*. [DOI](https://doi.org/10.1007/s10444-021-09843-0)
