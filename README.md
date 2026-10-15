# Skeletonized Wave-Equation Dispersion Spectrum Inversion (SWD)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Description
Skeletonized Wave-Equation Dispersion Spectrum Inversion (SWD), is a cutting-edge method in geophysical research aimed at obtaining a robust and reliable near-surface S-wave velocity structure. This method leverages a skeletal inversion framework that avoids traditional full waveform inversion's susceptibility to cycle-skipping by implementing a smooth gradient approximation between the dispersion spectrum and the misfit function. This is achieved through the SoftMax approximation.

The technique innovatively derives the gradient of the misfit function with respect to the velocity model utilizing the chain rule and adjoint state method. This integration allows SWD to couple with the wave equation, enabling precise and stable S-wave velocity inversions. Unlike conventional methods, SWD does not depend on a layered assumption, thus enhancing lateral resolution significantly.

SWD capitalizes on the concept of skeletonizing complex surface wave arrivals into simpler forms—specifically, picked dispersion curves in the phase-velocity and frequency domains, akin to wave-equation traveltime tomography. These dispersion curves are primarily obtained from Rayleigh waves captured by vertical-component geophones. The misfit function itself is defined as the sum of the squared differences between the wavenumbers of the predicted and observed dispersion curves, reflecting the method's refined approach to accurately capturing subsurface velocity structures.

Link to this repository: https://github.com/zhangchang271/SWD

## Key Features
- Mid-high resolution compared to FWI.
- No need to estiamte source wavelet.
- The finite difference forward and inverse kernels are accelerated by C++.
- Elastic waves with topography and fractional source/receiver positions.
- Reverse-time source-wavefield reconstruction during adjoint propagation.

   
## Usage
Run `SWD.m` directly for the model test. It is the single model-test entry point and uses the source and receiver coordinates together with `dx` and `offset`; it does not convert the offset window through `M`, `m`, or shot-index boundary formulas.

Set `terrain` in the acquisition block: one surface z coordinate per model column at `x=(0:nx-1)*dx`, in the same physical units as `dx`, positive downward. `terrain=zeros(1,nx)` gives the original flat surface. Model rows follow the surface and stretch to the fixed bottom at `z=(nz-1)*dx`; model columns follow the horizontal grid deformation.

Continue supplying positions through `sx/sz` and `gx/gz`, including fractional coordinates, irregular spacing, and receiver depths. Their original convention is retained: physical `(x,z)` becomes `(1+x/dx,1+z/dx)`. The default geometry places sources and receivers on `terrain`; replace these arrays directly when using surveyed positions. The grid aligns separated receivers with nearby nodes and interpolates other positions, preserving gather order. Distances used by WD/SWD remain the original horizontal offsets.

The deformed grid uses the existing default `fd_order=22`, `source_type='w'`, and `isfs=1`. Its time step is adjusted internally for the grid deformation. The original uniform flat-grid stencils remain available at orders 22, 24, 26, and 28.

`SWD.m` prepares the fixed grid, receiver positions, and propagation settings once in `grid`. Each shot then uses:

```matlab
[u,w,forward] = staggerfd_eigen(vp,vs,s,sx,sz,grid);
[g_lambda,g_mu,illum] = e2drtm_eigen(single(residual),forward);
```

The output order remains horizontal data, vertical data, and the information needed by the corresponding adjoint call. A forward-only call uses `[u,w]`. The source length determines the number of time steps; model and grid settings are passed internally from forward to adjoint propagation.

Wavefield reconstruction is automatic. Forward propagation saves the final state and a thin boundary strip, then reconstructs the source field alongside the adjoint. Boundary records use an internal temporary file, removed after the corresponding adjoint call; full volume snapshots are no longer stored through time. `nt_wf` retains its role as the gradient accumulation interval. No reconstruction switch is needed.

`SWD.m` builds its two active MEX kernels before starting parallel workers. A configured MATLAB C++ compiler and Eigen 3 headers are required (`/usr/include/eigen3`, or `core/eigen3`). These kernels are built locally because their interface now includes the fixed grid and reconstruction data.

The former index-window main program, single-precision example, and live script are available on the [`legacy-master` branch](https://github.com/zhangchang271/SWD/tree/legacy-master). `SWDlive.pdf` documents the intermediate steps of the original program.

The fieldexamples folder contains two field examples, including processed dispersion curves, the main run file, and an mlx file that explains parameter selection in detail.


The distance workflow switches between WD and SWD by changing one call in `SWD.m`:
- **WD method (default)**: call `weight_dataAD`, which uses physical source-receiver distances.
- **SWD method**: uncomment the `ADWDgrad_w` call. Its Radon helpers return logical receiver masks in the original gather order.

The distance-aware Radon functions are `core/RTrADx.m` and `core/RTlADx.m`. Receiver selection is shared by these functions and `core/weight_dataAD.m` through `core/offset_traces.m`.

Run `addpath('tests'); test_fd` from the repository root to check the original flat-grid reference, irregular receiver sampling, a mapped-grid model derivative, reconstruction against full wavefield storage, and analytical Rayleigh-wave speeds on horizontal and inclined surfaces. Test compilation and full trajectories are confined to temporary storage.

## Result
![fig1.png](fig1.png)

## License
SWD is distributed under the GNU General Public License v3.0. See the `LICENSE` file for more details.

## Contact
This program was written by Zhang Chang under the supervision of Professor Li Jing from Jilin University. If you have any questions, please contact:
- **Zhang Chang**: zhang.chang271@gmail.com
[![Email](https://img.shields.io/badge/Email-zhang.chang271@gmail.com-blue)](mailto:zhang.chang271@gmail.com) Wechat: zc13604682616
- **Li Jing**: inter.lijing@gmail.com
  [![Email](https://img.shields.io/badge/Email-inter.lijing@gmail.com-red)](mailto:inter.lijing@gmail.com) 

## Citation
If you use this package in your own research, please cite the following papers:

Chang Zhang, Jing Li, Sherif Hanafy, Lige Bai, Hui Liu, Huaqing Cao, Yibo Wang; Wave-equation skeletonized inversion of dispersion curves and spectra. Geophysics 2026; doi: https://doi.org/10.1190/GEO-2024-0762
