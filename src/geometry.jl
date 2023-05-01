
try
    using Gmsh: gmsh
catch
    using gmsh
end

include("definitions/constants.jl")


gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)


function gmsh_add_rectangle(mid, width, height, lc)
    geo = gmsh.model.geo;
    
    # Corner points
    p1 = geo.addPoint(mid[1] - width / 2, mid[2] - height / 2, 0, lc);
    p2 = geo.addPoint(mid[1] + width / 2, mid[2] - height / 2, 0, lc);
    p3 = geo.addPoint(mid[1] + width / 2, mid[2] + height / 2, 0, lc);
    p4 = geo.addPoint(mid[1] - width / 2, mid[2] + height / 2, 0, lc);
    points = [p1, p2, p3, p4];
    
    # Lines
    l1 = geo.addLine(p1, p2);
    l2 = geo.addLine(p2, p3);
    l3 = geo.addLine(p3, p4);
    l4 = geo.addLine(p4, p1);
    lines = [l1, l2, l3, l4];
    
    # Curve loop
    loop = geo.addCurveLoop(lines);
    
    return loop, lines, points;
end



function gmsh_add_rectangle(mid, width, height, lc, radius)
    geo = gmsh.model.geo;
    
    # Corner points
    p1 = geo.addPoint(mid[1] - width / 2 + radius, mid[2] - height / 2, 0, lc);
    p2 = geo.addPoint(mid[1] + width / 2 - radius, mid[2] - height / 2, 0, lc);
    p3 = geo.addPoint(mid[1] + width / 2         , mid[2] - height / 2 + radius, 0, lc);
    p4 = geo.addPoint(mid[1] + width / 2         , mid[2] + height / 2 - radius, 0, lc);
    p5 = geo.addPoint(mid[1] + width / 2 - radius, mid[2] + height / 2, 0, lc);
    p6 = geo.addPoint(mid[1] - width / 2 + radius, mid[2] + height / 2, 0, lc);
    p7 = geo.addPoint(mid[1] - width / 2         , mid[2] + height / 2 - radius, 0, lc);
    p8 = geo.addPoint(mid[1] - width / 2         , mid[2] - height / 2 + radius, 0, lc);
    
    c1 = geo.addPoint(mid[1] - width / 2 + radius, mid[2] - height / 2 + radius, 0, 1);
    c2 = geo.addPoint(mid[1] + width / 2 - radius, mid[2] - height / 2 + radius, 0, 1);
    c3 = geo.addPoint(mid[1] + width / 2 - radius, mid[2] + height / 2 - radius, 0, 1);
    c4 = geo.addPoint(mid[1] - width / 2 + radius, mid[2] + height / 2 - radius, 0, 1);
    points = [p1, p2, p3, p4, p5, p6, p7, p8];
    
    # Lines
    l1 = geo.addLine(p1, p2)
    l2 = geo.addCircleArc(p2, c2, p3)
    l3 = geo.addLine(p3, p4)
    l4 = geo.addCircleArc(p4, c3, p5)
    l5 = geo.addLine(p5, p6)
    l6 = geo.addCircleArc(p6, c4, p7)
    l7 = geo.addLine(p7, p8)
    l8 = geo.addCircleArc(p8, c1, p1)
    lines = [l1, l2, l3, l4, l5, l6, l7, l8]
    
    # Curve loop
    loop = geo.addCurveLoop(lines);
    
    return loop, lines, points;
end


gmsh.model.add("transformer")
print("Added transformer.\n")

geo = gmsh.model.geo;

## Enclosure
enclosure_lp, enclosure_lines, _ = gmsh_add_rectangle([0, 0.5 * (hencl1 + hencl2)], wencl, hencl1 - hencl2, lc1)
print("Done with enclosure.\n")

## Core 
core_lp, _, _  = gmsh_add_rectangle([0, 0], wcore, hcore, lc1)     # Core outline
cgap1_lp, _, _ = gmsh_add_rectangle([-mgap, 0], wgap, hgap, lc2, 1.5e-2)   # Gap left
cgap2_lp, _, _ = gmsh_add_rectangle([+mgap, 0], wgap, hgap, lc2, 1.5e-2)   # Gap right
print("Done with core.\n")

## HV windings
xm = mgap + wgap / 2 + (wcore / 2 - mgap - wgap / 2) / 2;

hv1l_lp, _, _ = gmsh_add_rectangle([-xm - mwhv, 0], wwhv, hwhv, lc3)
hv1r_lp, _, _ = gmsh_add_rectangle([-xm + mwhv, 0], wwhv, hwhv, lc3)
hv2l_lp, _, _ = gmsh_add_rectangle([    - mwhv, 0], wwhv, hwhv, lc3)
hv2r_lp, _, _ = gmsh_add_rectangle([    + mwhv, 0], wwhv, hwhv, lc3)
hv3l_lp, _, _ = gmsh_add_rectangle([+xm - mwhv, 0], wwhv, hwhv, lc3)
hv3r_lp, _, _ = gmsh_add_rectangle([+xm + mwhv, 0], wwhv, hwhv, lc3)
print("Done with HV windings.\n")

## LV windings
xm = mgap + wgap / 2 + (wcore / 2 - mgap - wgap / 2) / 2;

lv1l_lp, _, _ = gmsh_add_rectangle([-xm - mwlv, 0], wwlv, hwlv, lc4)
lv1r_lp, _, _ = gmsh_add_rectangle([-xm + mwlv, 0], wwlv, hwlv, lc4)
lv2l_lp, _, _ = gmsh_add_rectangle([    - mwlv, 0], wwlv, hwlv, lc4)
lv2r_lp, _, _ = gmsh_add_rectangle([    + mwlv, 0], wwlv, hwlv, lc4)
lv3l_lp, _, _ = gmsh_add_rectangle([+xm - mwlv, 0], wwlv, hwlv, lc4)
lv3r_lp, _, _ = gmsh_add_rectangle([+xm + mwlv, 0], wwlv, hwlv, lc4)
print("Done with LV windings.\n")


## Surfaces
geo.addPlaneSurface([enclosure_lp, core_lp, hv1l_lp, lv1l_lp, hv3r_lp, lv3r_lp], 1)
geo.addPlaneSurface([core_lp, cgap1_lp, cgap2_lp], 2)
geo.addPlaneSurface([cgap1_lp, hv1r_lp, lv1r_lp, hv2l_lp, lv2l_lp], 3)
geo.addPlaneSurface([cgap2_lp, hv2r_lp, lv2r_lp, hv3l_lp, lv3l_lp], 4)

geo.addPlaneSurface([hv1l_lp], 5)
geo.addPlaneSurface([hv1r_lp], 6)
geo.addPlaneSurface([hv2l_lp], 7)
geo.addPlaneSurface([hv2r_lp], 8)
geo.addPlaneSurface([hv3l_lp], 9)
geo.addPlaneSurface([hv3r_lp], 10)

geo.addPlaneSurface([lv1l_lp], 11)
geo.addPlaneSurface([lv1r_lp], 12)
geo.addPlaneSurface([lv2l_lp], 13)
geo.addPlaneSurface([lv2r_lp], 14)
geo.addPlaneSurface([lv3l_lp], 15)
geo.addPlaneSurface([lv3r_lp], 16)
print("Done with surfaces.\n")


geo.synchronize()
print("Synchronizing.\n")


## Physical Groups
geo.addPhysicalGroup(2, [1, 3, 4], 1)   # Transformer oil
geo.addPhysicalGroup(2, [2], 2)         # Core
geo.addPhysicalGroup(2, [5], 3)         # HV winding phase 1 left
geo.addPhysicalGroup(2, [6], 4)         # HV winding phase 1 right
geo.addPhysicalGroup(2, [7], 5)         # HV winding phase 2 left
geo.addPhysicalGroup(2, [8], 6)         # HV winding phase 2 right
geo.addPhysicalGroup(2, [9], 7)         # HV winding phase 3 left
geo.addPhysicalGroup(2, [10], 8)        # HV winding phase 3 right

geo.addPhysicalGroup(2, [11], 9)        # LV winding phase 1 left
geo.addPhysicalGroup(2, [12], 10)       # LV winding phase 1 right
geo.addPhysicalGroup(2, [13], 11)       # LV winding phase 2 left
geo.addPhysicalGroup(2, [14], 12)       # LV winding phase 2 right
geo.addPhysicalGroup(2, [15], 13)       # LV winding phase 3 left
geo.addPhysicalGroup(2, [16], 14)       # LV winding phase 3 right

geo.addPhysicalGroup(2, [5, 6, 7, 8, 9, 10], 15)      # HV windings
geo.addPhysicalGroup(2, [11, 12, 13, 14, 15, 16], 16) # LV windings

geo.addPhysicalGroup(1, enclosure_lines, 1)  # Enclosure boundary

gmsh.model.setPhysicalName(2, 1, "Oil")
gmsh.model.setPhysicalName(2, 2, "Core")
gmsh.model.setPhysicalName(2, 3, "HV1l")
gmsh.model.setPhysicalName(2, 4, "HV1r")
gmsh.model.setPhysicalName(2, 5, "HV2l")
gmsh.model.setPhysicalName(2, 6, "HV2r")
gmsh.model.setPhysicalName(2, 7, "HV3l")
gmsh.model.setPhysicalName(2, 8, "HV3r")
gmsh.model.setPhysicalName(2, 9, "LV1l")
gmsh.model.setPhysicalName(2, 10, "LV1r")
gmsh.model.setPhysicalName(2, 11, "LV2l")
gmsh.model.setPhysicalName(2, 12, "LV2r")
gmsh.model.setPhysicalName(2, 13, "LV3l")
gmsh.model.setPhysicalName(2, 14, "LV3r")
gmsh.model.setPhysicalName(2, 15, "HV windings")
gmsh.model.setPhysicalName(2, 16, "LV windings")

gmsh.model.setPhysicalName(1, 1, "Enclosure")
print("Done with physical groups.\n")


# Generate Mesh
geo.synchronize()
gmsh.model.mesh.generate(2)
print("Done with generating mesh.\n")

try
    gmsh.write("geo/transformer_stedin.msh")
    print("Mesh was successfully written.\n")

catch
    print("Error!")
end