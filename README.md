# Exact and Efficient Polyhedral Envelope Containment Check

![Build](https://github.com/wangbolun300/fast-envelope/workflows/Build/badge.svg)

![](bunny.jpg)
This is the implementation of the paper "Exact and Efficient Polyhedral Envelope Containment Check". Our algorithm conservatively tells you if the distance between a query triangle and a triangle mesh is within a user-specific parameter $\epsilon$, and thus can be used to conservatively/safely reject a triangle that is too far from a triangle mesh. 

Using a polyhedral envelope presentation, our distance bound is a little bit tighter than the traditional Hausdorff distance bound, and our polyhedral-envelope check is exact! It means that if our algorithm returns "true" when checking if a triangle is within $\epsilon$ distance from a triangle mesh, it is also within a Hausdorff distance of $\epsilon$, without any uncertainty caused by numerical instability. Our method is significantly faster than the traditional sampling-based Hausdorff distance prediction methods when $\epsilon$ is small, e.g., $10^{-4}$ of the diagonal length of the axis-aligned bounding box of the reference triangle mesh.

If you use our code, please cite our paper
```bibtex
@article{Wang:2020:FE,
    title={Exact and Efficient Polyhedral Envelope Containment Check},
    author={Bolun Wang and Teseo Schneider and Yixin Hu and Marco Attene and Daniele Panozzo},
    journal = {ACM Trans. Graph.},
     volume = {39},
     number = {4},
     month = jul,
     year = {2020},
     publisher = {ACM}
}
```
Please click [HERE](https://cims.nyu.edu/gcl/papers/2020-Fast-Envelope.pdf) to download the paper.
This is the link to our talk on SIGGRAPH 2020 [https://www.youtube.com/watch?v=_Vm61nlxyBI](https://www.youtube.com/watch?v=_Vm61nlxyBI).

**Implementation in CGAL**

A partial reimplementation of the algorithm in this repository, which does not use the indirect predicates, is available in CGAL 5.3 (https://doc.cgal.org/latest/Polygon_mesh_processing/index.html#title36). We still recommend you use our implementation, since our code has been fully tested and evaluated in multiple remeshing applications, including [fTetWild](https://github.com/wildmeshing/fTetWild).

## Important Note
There is a compiler flag that is required to ensure the correctness of the algorithm.
The flag is not available on Clang. The code has been tested on GCC and the Windows compiler.


# Installation via CMake
 - Clone our repository in your dependency folder (or add it as a submodule)
 - Add this in your main `CMakeLists.txt` file, `add_subdirectory` pointing to the directory where you cloned this repository
 - link your target with our library `target_link_libraries(<your-target> PUBLIC FastEnvelope)`

 ## Note
 Our library requires standard predicates to work; by default, we use the fast predicates inside [Geogram](http://alice.loria.fr/software/geogram/doc/html/index.html). If you want to avoid having Geogram as a dependency, you can disable it by setting `FAST_ENVELOPE_WITH_GEOGRAM_PSM_PREDICATES` to `ON`. The code will be slower.

 # Usage
  - Include `#include <fastenvelope/FastEnvelope.h>`
  - Initialize the envelope checker `FastEnvelope(const std::vector<Vector3>& m_ver, const std::vector<Vector3i>& m_faces, const Scalar eps);` with vertices, connectivity, and envelope size.
  - Call one of the `is_outside` functions with a triangle, point, or segment.


 # Testing
 We also provide an executable target `FastEnvelope_bin` that can be used for benchmarking

 You can run it by:
```bash
./FastEnvelope_bin ./queries/<INPUT>_envelope_log.csv ./ftetwild_queries/<INPUT> <OUTPUT> 1e-3 1 ours.
```

# Data
All data used in our paper can be downloaded from [https://archive.nyu.edu/handle/2451/61221](https://archive.nyu.edu/handle/2451/61221).
