using gmsh

using LinearAlgebra, SparseArrays
using WriteVTK

using BenchmarkTools

include("constants.jl")
include("load_geometry.jl")


const MESH_LOCATION = "geo/transformer_stedin.msh"



function fem2d(mshdata, sourceperelement, reluctivityperelement, conductivityperelement, omega)
    A = spzeros(Complex{Float64}, mshdata.nnodes, mshdata.nnodes)
    f = zeros(Complex{Float64}, mshdata.nnodes, 1)

    xnode = mshdata.xnode;
    ynode = mshdata.ynode;
    
    # Perform a loop over the elements
    for (element_id, nodes) in enumerate(mshdata.elements)
        # Retrieve global numbering of the local nodes of the current element
        node1_id = nodes[1]; node2_id = nodes[2]; node3_id = nodes[3];

        # Retrieve the x and y coordinates of the local nodes of the current element
        xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
        ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];

        # Compute surface area of the current element
        x12 = xnode2 - xnode1; x13 = xnode3-xnode1;
        y12 = ynode2 - ynode1; y13 = ynode3-ynode1;
        area_id = x12*y13 - x13*y12; area_id = abs(area_id)/2

        # Compute local vector contribution floc of the current element
        floc = area_id/3*sourceperelement[element_id] * [1; 1; 1]

        # Compute local matrix contribution Aloc of the current element
        Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]] \ UniformScaling(1.);
        Emat[3,:] .= 0;
        Bloc = area_id*reluctivityperelement[element_id]*(transpose(Emat)*Emat);
        Cloc = 1im * area_id / 3 * conductivityperelement[element_id] * omega * Diagonal(ones(3));

        # Add local contribution to f and A
        f[nodes]       += floc;
        A[nodes,nodes] += (Bloc + Cloc);
    end

    # Handle the boundary conditions
    bnd_node_ids, _ = gmsh.model.mesh.getNodesForPhysicalGroup(1, 1);
    A[bnd_node_ids,:] .= 0;
    A[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)))
    f[bnd_node_ids] .= 0;
    
    # Compute the numerical solution
    u = A \ f
    
    return u;
end



function process(mshdata, u, sourceperelement, reluctivityperelement, conductivityperelement, omega)
    Bx = zeros(Complex{Float64}, mshdata.nelements);
    By = zeros(Complex{Float64}, mshdata.nelements);
    
    Jel = zeros(mshdata.nelements);
    
    # Perform a loop over the elements
    for (element_id, nodes) in enumerate(mshdata.elements)
        # Retrieve global numbering of the local nodes of the current element
        node1_id = nodes[1]; node2_id = nodes[2]; node3_id = nodes[3];
        
        # Get x and y coordinates of the three nodes
        xnode = mshdata.xnode; ynode = mshdata.ynode;
        xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
        ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];

        # Compute surface area of the current element
        x12 = xnode2 - xnode1; x13 = xnode3-xnode1;
        y12 = ynode2 - ynode1; y13 = ynode3-ynode1;
        area_id = x12*y13 - x13*y12; area_id = abs(area_id)/2

        # Calculate shape function parameters
        Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]] \ UniformScaling(1.);
    
        # Calculate Bx and By from the solution coefficients and the shape function parameters
        c = u[[node1_id, node2_id, node3_id]];
        Bx[element_id] = sum(c .* Emat[2,:]);
        By[element_id] = -sum(c .* Emat[1,:]);
        
        # Calculate eddy current loss
        sigma = conductivityperelement[element_id];
        Jel[element_id] = norm(sourceperelement[element_id] + omega * sigma * 1/3 * sum(c));
    end
    
    # H is related to B through the reluctivity
    Hx = reluctivityperelement' .* Bx;
    Hy = reluctivityperelement' .* By;
    
    # Energy is 0.5 * dot(B, H)
    Wm = 0.5 * (Bx .* Hx .+ By .* Hy);
    
    return (Bx,By), (Hx, Hy), Wm, Jel;
end




############################################################

# Linear FEM without eddy currents

############################################################


gmsh.open(MESH_LOCATION)
mshdata = get_mesh_data();



Ip = 0;       # Primary peak phase current
Is = 777.62;  # Secondary peak phase current
Np = 266;
Ns = 6;

omega = 2*pi*50;  # Frequency

# Calculate current density in the windings
Jp = Np * Ip / Awhv;
Js = Ns * Is / Awlv;

# Source current density J
#  One term for each of the windings, with a positive and negative part
#  Note the phase shift between the phases
sourcefunction(group_id) = Jp * exp(1im * 2pi/3) * (-1 * (group_id==3) + 1 * (group_id==4)) + 
                           Jp * (-1 * (group_id==5) + 1 * (group_id==6)) + 
                           Jp * exp(-1im * 2pi/3) * (-1 * (group_id==7) + 1 * (group_id==8)) + 
                           Js * exp(1im * 2pi/3) * (1 * (group_id==9) - 1 * (group_id==10)) +
                           Js * (1 * (group_id==11) - 1 * (group_id==12)) + 
                           Js * exp(-1im * 2pi/3) * (1 * (group_id==13) - 1 * (group_id==14));
sourceperelement = map(sourcefunction, mshdata.e_group);

# Relative permeability model
mu0 = 4e-7 * pi;
mur = 1000;       # Relative permeability of the core
reluctivityfunction(group_id) = (1 / mu0) + (1/(mu0*mur) - 1/mu0) * (group_id == 2)
reluctivityperelement = map(reluctivityfunction, mshdata.e_group);

# Conductivity
conductivityfunction(group_id) = 0;
conductivityperelement = map(conductivityfunction, mshdata.e_group);



# Calculate the vector potential
u = fem2d(mshdata, sourceperelement, reluctivityperelement, conductivityperelement, omega);

# Post-process for magnetic field and current density
B, H, Wm, Jel = process(mshdata, u, sourceperelement, reluctivityperelement, conductivityperelement, omega);
Bnorm = norm.(sqrt.(B[1].^2 + B[2].^2));


# Define nodes (points) and elements (cells)
points = [mshdata.xnode mshdata.ynode]';
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el) for el in mshdata.elements];

# Create VTK file structure using nodes and elements
vtkfile = vtk_grid("images/transformer/transformer1.vtu", points, cells);

# Store data in the VTK file
vtkfile["Az", VTKPointData()]   = norm.(u);
vtkfile["imA", VTKPointData()]  = imag.(u);
vtkfile["Bnorm", VTKCellData()] = Bnorm;
vtkfile["Jel", VTKCellData()]   = Jel;

# Save the file
outfiles = vtk_save(vtkfile);

