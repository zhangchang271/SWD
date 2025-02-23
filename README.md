# Skeletonized Wave-Equation Dispersion Spectrum Inversion (SWD)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Description
Skeletonized Wave-Equation Dispersion Spectrum Inversion (SWD), is a cutting-edge method in geophysical research aimed at obtaining a robust and reliable near-surface S-wave velocity structure. This method leverages a skeletal inversion framework that eschews traditional full waveform inversion's susceptibility to cycle-skipping by implementing a smooth gradient approximation between the dispersion spectrum and the misfit function. This is achieved through the SoftMax approximation.

The technique innovatively derives the gradient of the misfit function with respect to the velocity model utilizing the chain rule and adjoint state method. This integration allows SWD to couple with the wave equation, enabling precise and stable S-wave velocity inversions. Unlike conventional methods, SWD does not depend on a layered assumption, thus enhancing lateral resolution significantly.

SWD capitalizes on the concept of skeletonizing complex surface wave arrivals into simpler forms—specifically, picked dispersion curves in the phase-velocity and frequency domains, akin to wave-equation traveltime tomography. These dispersion curves are primarily obtained from Rayleigh waves captured by vertical-component geophones. The misfit function itself is defined as the sum of the squared differences between the wavenumbers of the predicted and observed dispersion curves, reflecting the method's refined approach to accurately capturing subsurface velocity structures.

## Key Features
- Mid-high resolution compared to FWI.
- No need to estiamte source wavelet.
- Accelerated by C++: The finite difference forward and inverse kernels used in SWD are accelerated using C++.
- Suitable for elastic wave and flat surface conditions.

   
## Usage
run SWD.m

## Result
![fig1.png](fig1.png)

## License
SWD is distributed under the GNU General Public License v3.0. See the `LICENSE` file for more details.

## Contact
- **Zhang Chang**: zhangchang23@mails.jlu.edu.cn
[![Email](https://img.shields.io/badge/Email-zhangchang23@mails.jlu.edu.cn-blue)](mailto:zhangchang23@mails.jlu.edu.cn)  Jilin University, China
- **Li Jing**: inter.lijing@gmail.com
  [![Email](https://img.shields.io/badge/Email-inter.lijing@gmail.com-red)](mailto:inter.lijing@gmail.com)   Jilin University, China

## Citation
If you use this package in your own research, please cite the following papers:

- Zhang, C., J. Li, and H. Cao, 2024, Skeletonized wave-equation dispersion spectrum inversion: GEOPHYSICS (under review)
