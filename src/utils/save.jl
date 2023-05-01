
"""Save result"""
function save(file_name, mshdata, u, B, H, Wm, Jel)

    # Define nodes (points) and elements (cells)
    points = [mshdata.xnode mshdata.ynode]';
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el) for el in mshdata.elements];

    # Create VTK file structure using nodes and elements
    vtkfile = vtk_grid(string(OUTPUT_LOCATION, file_name), points, cells);

    # Store data in the VTK file
    vtkfile["Az", VTKPointData()]   = norm.(u);
    vtkfile["imA", VTKPointData()]  = imag.(u);
    vtkfile["Bnorm", VTKCellData()] = norm.(sqrt.(B[1].^2 + B[2].^2));
    vtkfile["B_vec", VTKCellData()] = real.(B)
    vtkfile["Jel", VTKCellData()]   = Jel;

    # Save the file
    print("Saving result in a file...")
    outfiles = vtk_save(vtkfile);
    println(" Done.")

end