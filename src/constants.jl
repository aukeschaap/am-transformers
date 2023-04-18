

# Physical constants
# Enclosure dimensions (enclosure is not centered, the core is)
hencl1 = 65.5e-2;   # Height above the x-axis
hencl2 = -53.5e-2;  # Height below the y-axis
wencl  = 104e-2;    # Width

# Core dimensions
wcore = 84e-2;
hcore = 100e-2;

# Core gap dimensions (left and right are identical)
wgap = 17e-2;
hgap = 76e-2;
mgap = 17e-2;

# HV winding dimensions (all phases left/right are identical)
wwhv = 3e-2;
hwhv = 74e-2;
mwhv = 14.75e-2;
Awhv = wwhv * hwhv;

# LV winding dimensions (all phases left/right are identical)
wwlv = 2e-2;
hwlv = 74e-2;
mwlv = 11.25e-2;
Awlv = wwlv * hwlv;

# Mesh densities
lc1 = 2e-2;      # Enclosure & core outer
lc2 = 1e-2;      # Core inner
lc3 = 1e-2;      # HV windings
lc4 = 1e-2;      # LV windings