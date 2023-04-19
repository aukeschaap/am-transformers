# Modeling magnetic field in a power transformer
By Philip Soliman (4945255) and Auke Schaap (4457919).

This is a project for the course WI4204 Advanced Modeling at the TU Delft, in combination with STEDIN Rotterdam. The project is supervised by Domenico Lahaye.


> ⚠️ Make sure that `gmsh` is setup correctly if you are using Windows! Follow [this](https://github.com/ziolai/finite_element_electrical_engineering/blob/main/extended-lab-sessions/gmsh/Mesh-Generation-using-Gmsh.ipynb) guide by dr. Domenico Lahaye. The alternative is to use [WSL](https://learn.microsoft.com/en-us/windows/wsl/about).


## Structure

The project is separated in the following way:
- `doc` contains relevant (pdf) documents, e.g. Max's Thesis
- `mesh` contains all the (generated) meshes
- `report` has the LaTeX files that are used to build the report
- `src` holds the Julia source code


# Useful links

#### Project repo
https://github.com/ziolai/finite_element_electrical_engineering (project description)

https://github.com/ziolai/finite_element_electrical_engineering/blob/main/lab-sessions/6-lab-session.ipynb (time integration, section 4)

#### Repos of Gijs Lagerweij
https://github.com/gijswl/ee4375_fem (general FEM)

https://github.com/gijswl/ee4375_fem_ta (distribution transformer)

#### Gmsh repo
https://github.com/JuliaFEM/Gmsh.jl (gmsh for julia)

https://github.com/ziolai/finite_element_electrical_engineering/blob/main/extended-lab-sessions/gmsh/Mesh-Generation-using-Gmsh.ipynb (gmsh for julia workaround + intro to gmsh)

#### Gmsh documentation 
https://gmsh.info/doc/texinfo/gmsh.html#index-Command_002dline-options

#### DifferentialEquations julia library
https://docs.sciml.ai/DiffEqDocs/stable/examples/diffusion_implicit_heat_equation/ (time integration of diffusion equation using ODEproblem method)