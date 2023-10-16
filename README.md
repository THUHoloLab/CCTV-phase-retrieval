# Compressive phase retrieval via constrained complex total variation regularization (CCTV)
**Authors:** [Yunhui Gao](https://github.com/Yunhui-Gao) (gyh21@mails.tsinghua.edu.cn) and Liangcai Cao (clc@tsinghua.edu.cn)

<p align="center">
<img src="docs/diagram.png">
</p>

<p align="center"> <strong>Figure 1</strong>. (a) Schematic of the in-line holographic imaging system. (b) Captured raw hologram of a transparent Fresnel zone plate. Scale bar 1 mm. (c) Retrieved phase distribution. (d) Rendered surface height profile.</p>

## Requirements
Matlab 2019a or newer. Older visions may be sufficient but have not been tested.

## Quick Start
- **Phase retrieval using simulated data.** Run [`demo_sim.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_sim.m) with default parameters.
- **Phase retrieval using experimental data.** First follow the instruction [here](https://github.com/THUHoloLab/CCTV-phase-retrieval/tree/master/data/experiment) to download the data. Then run [`demo_exp.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_exp.m) with default parameters.
- **Try on your own experiment data.** Prepare a hologram and an optional reference image, run [`preprocessing.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/preprocessing.m) and set the experiment parameters (e.g. pixel size, wavelength, and sample-to-sensor distance). Then run [`demo_exp.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_exp.m) and see how it works.


## Accelerated implementations
The basic demo codes provide intuitive and proof-of-concept implementations for beginners, but are far from efficient. To facilitate faster reconstruction, we provide an optimized version based on CPU or GPU, which can be found at [`demo_sim_fast.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_sim_fast.m) and [`demo_exp_fast.m`](https://github.com/THUHoloLab/CCTV-phase-retrieval/blob/master/main/demo_exp_fast.m). Figure 2 and Table 1 show the runtime (200 iterations) for different image dimensions. The results are obtained using a laptop computer with Intel&reg; Core&trade; i7-12700H (2.30 GHz) CPU and Nvidia GeForce RTX&trade; 3060 GPU.

<p align="center">

|  Image dimension    | CPU runtime (s) | GPU runtime (s) |
|  :----:             | :----:          | :----:          |
|  128 $\times$ 128   | 0.673           | 0.704           |
|  256 $\times$ 256   | 2.76            | 0.824           |
|  512 $\times$ 512   | 8.76            | 1.25            |
|  1024 $\times$ 1024 | 31.8            | 3.67            |
|  2048 $\times$ 2048 | 130.8           | 13.2            |

</p>
<p align="center"> <strong>Table 1</strong>. Runtimes (for 200 iterations) using GPU and CPU for different image dimensions.</p>


<p align="center">
<img src="docs/runtime.png" style="height: 300px;">
</p>

<p align="center"> <strong>Figure 2</strong>. Runtimes (for 200 iterations) using GPU and CPU for different image dimensions.</p>

## Theories and References
For algorithm derivation and implementation details, please refer to our paper:

- [Yunhui Gao and Liangcai Cao, "Iterative projection meets sparsity regularization: towards practical single-shot quantitative phase imaging with in-line holography," Light: Advanced Manufacturing 4, 6 (2023).](https://www.light-am.com/article/doi/10.37188/lam.2023.006)
