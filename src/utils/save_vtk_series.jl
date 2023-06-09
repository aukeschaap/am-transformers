import JSON
using Formatting
include("save_vtk.jl")
include("solution.jl")

function save_vtk_series(name :: String, mesh_data:: MeshData, sol, source_per_element, reluctivity_per_element, conductivity_per_element)
        
        fspec = FormatSpec("<d") # flush left, decimal integer

        name *= "_" * fmt(fspec, 0)
        solution_folder = joinpath(OUTPUT_LOCATION, name)

        #check if folder exists
        passed = false
        i = 0
        while !passed
                i += 1
                if isdir(solution_folder)
                        exists = true 
                        name = replace(name, fmt(fspec, i-1) => fmt(fspec, i))
                        solution_folder = joinpath(OUTPUT_LOCATION, name)        
                        continue                
                end
                passed = true
        end

        mkdir(solution_folder)
        time_series_fn = joinpath(solution_folder, name * "_series")
        time_series = paraview_collection(time_series_fn)

        n_time_steps = length(sol.t)
        len = 0
        for i in 1:n_time_steps
                t = sol.t[i]
                u = sol.u[i]
                fspec = FormatSpec("<d") # flush left, decimal integer
                formatted_time = fmt(fspec, floor(t * 1e6)) #convert(Int64, t * 1e6)) <- causes inexactly representable error
                time_stamp = "t=" * formatted_time * "μs"
                file_name_t = name * "_" * time_stamp
                progress = "saved file: " * file_name_t * " total progress : " * string(round(i/n_time_steps*100, digits=2)) * "%\r"
                if length(progress) > len
                        len = length(progress)
                end
                print(progress)

                # save the vtk file
                B, H, Wm, Jel = solution(mesh_data, u, source_per_element, reluctivity_per_element, conductivity_per_element)
                vtk_path = joinpath(name, file_name_t * ".vtu")  # relative to OUTPUT_LOCATION
                vtk_file = save_vtk(vtk_path, mesh_data, u, B, H, Wm, Jel)

                # save the file name and time
                collection_add_timestep(time_series, vtk_file, t)
        end
        # clear progress bar
        for i in 1:len
                print(" ")
        end     
        print("\r")

        # save time series data   
        vtk_save(time_series)
end