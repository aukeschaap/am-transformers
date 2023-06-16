# GPU

We tried to use the GPU to calculate. The conclusion of the experiment is that it is not easy to use the GPU to speed up your code. You will either:
- Need to write your own kernel functions that run entirely on the GPU, and not transfer data between the GPU and the CPU.
- Need to use a package that can benefit from using the GPU, i.e. DifferentialEquations.jl.



### Links
https://www.intel.com/content/www/us/en/developer/articles/technical/vs-code-wsl2-and-oneapi-cross-platform-development.html (guide to install oneAPI)

https://juliagpu.org/

https://github.com/JuliaGPU/oneAPI.jl (oneAPI julia package)


https://github.com/intel/compute-runtime (required, including the right driver)

https://github.com/intel/compute-runtime/blob/master/WSL.md (info regarding wsl)

https://www.intel.com/content/www/us/en/download/19344/intel-graphics-windows-dch-drivers.html (example link for right drivers, actual link is located in wsl.md)

https://github.com/SciML/DiffEqGPU.jl (diff eq package for gpu)