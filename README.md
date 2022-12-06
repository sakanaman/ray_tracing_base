# Ray Tracing Base
Renderer for testing several techniques in Ray Tracing

![ibl_ridaisai](https://user-images.githubusercontent.com/42662735/84959115-b0760c00-b139-11ea-8b38-21186c92829d.png)

## Features
- Loose Coupling Libralies
    - Easily integrate them into your renderer
    - Easily test
- Cmake Build

## How To Build
### macOS 
```
mkdir build
cd build
cmake .. -Dembree_DIR=~/ray_tracing_base/ext/embree-3.13.2.x86_64.macosx/lib/cmake/embree-3.13.2
make
```
### Linux 
```
mkdir build
cd build
cmake .. -Dembree_DIR=~/ray_tracing_base/ext/embree-3.13.2.x86_64.linux/lib/cmake/embree-3.13.2
make
```

## Coming Soon ...
* Correspond to OpenVDB
* Implement Scene Graph
* More Light Transport Algorithms(BDPT, MCMC ...)
