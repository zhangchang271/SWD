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

Open **[SWD.m](SWD.m)** in MATLAB, edit the parameters at the top, and run it. This is the single entry point for the model example and contains the complete workflow.

Set `method = 'WD'` (default) or `method = 'SWD'` in the same script. WD uses `weight_dataAD`; SWD uses `ADWDgrad_w`.

`sx` and `gx` are grid coordinates; `dx` converts their difference to metres. `offset` sets the maximum source-receiver distance, and `minimum_traces` sets the minimum number of receivers on each side. Radon transforms and weighting share the selection rule in `core/offset_traces.m`.

The run leaves the final model `vs`, model history `vs_all`, gradients `dk_vs_all`, misfit history `residual`, and distance windows `offsets` in the MATLAB workspace.

`fieldexamples/` contains the two field datasets and their example workflows. `SWDlive.pdf` is an explanation of the original implementation. The original scripts and executable live tutorial are available on [legacy-master](https://github.com/zhangchang271/SWD/tree/legacy-master).

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
