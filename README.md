# ACP
(A)daptive (C)ontrolled (P)erturbation: a numerical library for implementing Computational Geometry algorithms.

# Requirements

All of the code in this repository has been tested under Ubuntu 18.04 with the
following requirements:

- General
  - cmake >= 3.11
  - pkg-config
  - pthread
  - gcc or clang with c++17 support

- Numerics
  - mpfr
  - gmp
  - lapack

- Visualization
  - OpenGL
  - GLFW3
  - GLUT
  - GLEW
  - GLM
  - VTK
  - freetype
  - z
  - png

To install all of the necessary packages at once in Ubuntu:
```
sudo apt install build-essential cmake pkg-config libopengl0 libglfw3-dev libglm-dev \
libglew-dev freeglut3-dev libgflags-dev libgflags2.2 libvtk7-dev libmpfr-dev libgmp-dev \
liblapack-dev liblapacke-dev libfreetype6-dev libpng-dev libpng16-16 libzip4
```

  # Building

  `cmake -S . -B build`

  `cmake --build build -j --target=surface_visualization_vtk`

  # Examples
  ## 4 Spheres
  `./build/bin/surface_visualization_vtk --surface=spheres --extent=2 --num_spheres=4`

  ## Torus and 2 spheres
  `./build/Debug/bin/surface_visualization_vtk --surface=torus,spheres --num_spheres=2 --mesh_epsilon=1e-2 --extent=2`

  ## Torus and Cayley
  `./build/Debug/bin/surface_visualization_vtk --surface=torus,cayley --mesh_epsilon=1e-2 --extent=2`

  ## Cayley and 2 spheres
  `./build/Debug/bin/surface_visualization_vtk --surface=cayley,spheres --mesh_epsilon=1e-2 --extent=2 --num_spheres=2`

  ## Torus
  `./build/bin/surface_visualization_vtk --surface=torus`

  ## Cayley
  `./build/Debug/bin/surface_visualization_vtk --surface=cayley --extent=3 --mesh_epsilon=1e-1`

  ## Clebsch
  `./build/Debug/bin/surface_visualization_vtk --surface=clebsch --extent=5 --angle=0.7 --mesh_epsilon=1e-1`

  ## Steiner
  `./build/Debug/bin/surface_visualization_vtk --surface=steiner --mesh_epsilon=1e-2`

  ## Hunt
  `./build/Debug/bin/surface_visualization_vtk --surface=hunt --extent=3 --mesh_epsilon=1e-2`