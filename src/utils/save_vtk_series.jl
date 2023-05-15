import JSON
using Formatting
include("save_vtk.jl")
include("solution.jl")

function save_vtk_series(name :: String, mesh_data:: MeshData, sol :: ODESolution, source_per_element, reluctivity_per_element, conductivity_per_element)
        
        n_time_steps = length(sol.t)

        vtk_series_dict = Dict(
                "file_series_version" => 1.0,
                "files" => Array{Dict}(undef, n_time_steps)
        )

        fspec = FormatSpec("<d") # flush left, decimal integer

        name *= "_" * fmt(fspec, 0)
        solution_folder = joinpath(OUTPUT_LOCATION, name)

        #check if folder exists
        passed = false
        i = 0
        while !passed && i < 10
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

        len = 0
        for i in 1:n_time_steps
                t = sol.t[i]
                u = sol[i]
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
                save_vtk(vtk_path, mesh_data, u, B, H, Wm, Jel)

                # save the file name and time in a dictionary
                vtk_series_dict["files"][i] = Dict(
                        "name" => file_name_t * ".vtu",
                        "time" => t
                );
        end
        # clear progress bar
        for i in 1:len
                print(" ")
        end     
        print("\r")

        # time series data as a json string      
        stringdata = JSON.json(vtk_series_dict)

        # write the json string to a *.vtk.series file
        vtk_series_fn = name * ".vtk.series"
        vtk_series_path = joinpath(solution_folder, vtk_series_fn)
        open(vtk_series_path, "w") do f
                write(f, stringdata)
        end
end