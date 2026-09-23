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
- Suitable for elastic wave and flat surface conditions.

   
## Usage
Run `SWD.m` for the model test. It launches `SWD_distance.m`, the current physical-distance workflow. `SWD_distance.m` uses the source and receiver coordinates together with `dx` and `offset`; it does not convert the offset window through `M`, `m`, or shot-index boundary formulas.

`SWD_legacy.m` preserves the former index-window implementation for comparison. `SWD_single.m` is the original single-precision example.

SWD.mlx is a MATLAB live script (similar to Jupyter) that details the intermediate steps of the SWD program. 

The fieldexamples folder contains two field examples, including processed dispersion curves, the main run file, and an mlx file that explains parameter selection in detail.


The distance workflow switches between WD and SWD by changing one call in `SWD_distance.m`:
- **WD method (default)**: call `weight_dataAD`, which uses physical source-receiver distances.
- **SWD method**: uncomment the `ADWDgrad_w` call. Its Radon helpers return logical receiver masks in the original gather order.

The distance-aware Radon functions are `core/RTrADx.m` and `core/RTlADx.m`. Receiver selection is shared by these functions and `core/weight_dataAD.m` through `core/offset_traces.m`.

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
