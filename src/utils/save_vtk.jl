
"""Save result"""
function save_vtk(file_name, mesh_data, u, B, H, Wm, Jel)

    # Define nodes (points) and elements (cells)
    points = [mesh_data.xnode mesh_data.ynode]';
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el) for el in mesh_data.elements];

    # Create VTK file structure using nodes and elements
    vtkfile = vtk_grid(string(OUTPUT_LOCATION, file_name), points, cells);

    # Store data in the VTK file
    vtkfile["Az", VTKPointData()]   = norm.(u);
    vtkfile["imA", VTKPointData()]  = imag.(u);
    vtkfile["Bnorm", VTKCellData()] = norm.(sqrt.(B[1].^2 + B[2].^2));
    vtkfile["B_vec", VTKCellData()] = real.(B)
    vtkfile["Jel", VTKCellData()]   = Jel;

    # Save the file
    outfiles = vtk_save(vtkfile);
    return vtkfile

end