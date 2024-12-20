# Curved Three-director Cosserat shells in *Stark*

This repository contains the source code used to produce the results for our paper:

> Löschner, F., Fernández-Fernández, J.A., Jeske, S.R. and Bender, J., 2024. **Curved Three-Director Cosserat Shells with Strong Coupling**. Computer Graphics Forum (CGF), DOI: 10.1111/cgf.15183.

Please also see the [entry on our website](https://animation.rwth-aachen.de/publication/0589/) for access to the paper, the supplemental document, video and more.
For newcomers to [**Stark**](https://github.com/InteractiveComputerGraphics/stark), please also check out the [upstream repository](https://github.com/InteractiveComputerGraphics/stark) without the paper specific additions.

<p align=center>
 <img src="docs/images/scroll.png" width="260">
  &nbsp;&nbsp;
 <img src="docs/images/dumpling.png" width="260">
  &nbsp;&nbsp;
 <img src="docs/images/strip.png" width="260">
</p>

## Overview

This repository is a fork of [**Stark**](https://github.com/InteractiveComputerGraphics/stark) [1], a C++ simulation platform that provides easy access to state-of-the-art methods to robustly solve simulations of rigid and deformable objects in a strongly coupled manner.
For most users who are *not* specifically interested in implementation details of the paper or reproducing experiments, we recommend to use the upstream [**Stark**](https://github.com/InteractiveComputerGraphics/stark) repository without the modifications of this project.

Notable additions on top of the main Stark repository to facilitate the experiments of the paper include:

 - The specific simulation "scenes" for all experiments of the paper which define material parameters, boundary conditions, scripting etc. See [`examples/main.cpp`](examples/main.cpp).
 - General FEM basis functions and quadrature rules used to discretize and evaluate the Cosserat models. In addition to the standard first- and second-order Lagrange triangle elements (`Tri3`, `Tri6`) presented in the paper, the code here is slightly more general and also contains third-order triangle elements (`Tri10`) as well as bilinear and biquadratic quadrilaterals (`Quad4`, `Quad9`). See [`fem_elements.h`](stark/src/models/fem_elements.h).
 - The implementation of the three-director Cosserat plate and shell models, previously presented by Nebel et al. 2023 [2], that we solve in an incremental potential context with strong coupling. See [`EnergyMicropolarShell.cpp`](stark/src/models/deformables/surface/EnergyMicropolarShell.cpp). 
 - An implementation of a Kirchhoff-Love shell model with support for arbitrary hyperelastic materials as introduced by Wen and Barbič 2023 [3] for qualitative comparisons. See [`EnergyTriangleStrainWen23.cpp`](stark/src/models/deformables/surface/EnergyTriangleStrainWen23.cpp).

## Build instructions

The project uses CMake as a build system and should work on Windows, Linux and macOS.
In-source builds are supported, but we recommend building in a subfolder:
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
cmake --build . --target examples --parallel 16
```
To run the experiments shown in the paper, execute the binary created in the `examples` folder:
```bash
cd ../
./examples/examples
```
By default, this runs the small metal strip twisting scene. For other scenes, uncomment the respective scene at the bottom of the [`main.cpp`](examples/main.cpp) file.

A few notes:
 - Some scenes use meshes found in the "`models`" folder of the repository. When running these scenes, the application assumes that the current working directory of the process is the root of this repository (i.e. the "`models`" folders is reachable as "`./models`").
 - On Linux and macOS, we use Clang for code generation at runtime (compilation of derivatives) which requires Clang to be installed and its executable to be in `PATH`.
   On Windows, MSVC is used for compilation at runtime which might require adapting the path in the `compiler_command` variable in the file [`Compilation.cpp`](stark/extern/symx/src/Compilation.cpp)
 - By default, this repository uses direct solvers as implemented by Eigen. Through the Eigen support modules it's possible to use external direct solvers for improved performance:
   - **Intel MKL:** By setting the CMake variable `STARK_ENABLE_MKL` to `ON` you can enable Intel MKL for the direct solvers (LU, LDLT). The MKL LU decomposition was also used for the timings and experiments in the paper.
     For the compilation with MKL to succeed, the appropriate environment variables have to be set.
     On Windows, this is usually achieved by starting the IDE or compilation from an "Intel oneAPI command prompt" that can be opened from a shortcut in the start menu after installation of MKL.
     On Linux, MKL should be detected automatically if it was installed from the official repositories. (Tested with MKL 2025.0.1 and some earlier versions)
   - **Apple Accelerate:** By setting the CMake variable `STARK_ENABLE_ACCELERATE` to `ON`, it is possible to use direct solvers from the Apple [Accelerate](https://developer.apple.com/documentation/accelerate/sparse_solvers?language=objc) framework on macOS provided through the [AccelerateSupport](https://eigen.tuxfamily.org/dox/group__AccelerateSupport__Module.html) module of Eigen. Currently, this is used only for the indefinite LDLT solver.
 - On Apple Silicon devices, it is required to manually set the CMake variable `STARK_ENABLE_AVX` to `OFF`.
 - The project uses OpenMP for parallelization. On macOS, you might need to install OpenMP with Homebrew and ensure that it can be found by CMake (e.g. by specifying `CMAKE_PREFIX_PATH=/opt/homebrew/opt/libomp`).

## References
 - [1] [Fernández-Fernández, J.A., Lange, R. and Laible, S. and Arras, K.O. and Bender, J., 2024. **STARK: A Unified Framework for Strongly Coupled Simulation of Rigid and Deformable Bodies with Frictional Contact**. 2024 IEEE International Conference on Robotics and Automation (ICRA).](https://ieeexplore.ieee.org/document/10610574)
 - [2] [Nebel, L.J., Sander, O., Bîrsan, M. and Neff, P, 2023. **A geometrically nonlinear Cosserat shell model for orientable and non-orientable surfaces: Discretization with geometric finite elements**. Computer Methods in Applied Mechanics and Engineering.](https://www.sciencedirect.com/science/article/pii/S0045782523004334)
 - [3] [Wen, J. and Barbič, J., 2023. **Kirchhoff-Love Shells with Arbitrary Hyperelastic Materials**. ACM SIGGRAPH Asia 2023.](https://dl.acm.org/doi/10.1145/3618405)

## License

The code in this repository (except for the files in the "`models`" folder that are explicitly mentioned below) is licensed under the Apache-2.0 license (in accordance with Stark itself at the time when this repository was forked). See [`LICENSE`](LICENSE).

The `armadillo_*` and `bunny_*` mesh variants in the "`models`" folder are variations of the [Armadillo and Stanford Bunny meshes, courtesy of the Stanford Computer Graphics Laboratory](http://graphics.stanford.edu/data/3Dscanrep/).

