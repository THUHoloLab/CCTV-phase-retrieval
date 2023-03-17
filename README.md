# Compressive phase retrieval via constrained complex total variation regularization (CCTV)
**Authors:** [Yunhui Gao](https://github.com/Yunhui-Gao) (gyh21@mails.tsinghua.edu.cn) and Liangcai Cao (clc@tsinghua.edu.cn)

<p align="center">
<img src="docs/diagram.png">
</p>

<p align="center"> Figure 1. (a) Schematic of the in-line holographic imaging system. (b) Captured raw hologram of a transparent Fresnel zone plate. Scale bar 1 mm. (c) Retrieved phase distribution. (d) Rendered surface height profile.</p>

## Requirements
Matlab 2019a or newer. Older visions may be sufficient but have not been tested.

## Usage
- **Phase retrieval using simulated data.** Run [`demo_sim.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_sim.m) with default parameters.
- **Phase retrieval using experimental data.** First follow the instruction [here](https://github.com/THUHoloLab/CCTV-phase-retrieval/tree/master/data/experiment) to download the data. Then run [`demo_exp.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_exp.m) with default parameters.
- **Try on your own experiment data.** Prepare a hologram and an optional reference image, run [`preprocessing.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/preprocessing.m) and set the experiment parameters (e.g. pixel size, wavelength, and sample-to-sensor distance). Then run [`demo_exp.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_exp.m) and see how it works.


## Theories and References
For algorithm derivation and implementation details, please refer to our paper:

- [Yunhui Gao and Liangcai Cao, "Iterative projection meets sparsity regularization: towards practical single-shot quantitative phase imaging with in-line holography," Light: Advanced Manufacturing 4, 6 (2023).](https://www.light-am.com/article/doi/10.37188/lam.2023.006)