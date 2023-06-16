# GPU

We tried to use the GPU to calculate. The conclusion of the experiment is that it is not easy to use the GPU to speed up your code. You will either:
- Need to write your own kernel functions that run entirely on the GPU, and not transfer data between the GPU and the CPU.
- Need to use a package that can benefit from using the GPU, i.e. DifferentialEquations.jl.