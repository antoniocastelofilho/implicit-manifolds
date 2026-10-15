import numpy as np
from fch_wrapper import FCH_CPU_Runner
import pyvista as pv

runner = FCH_CPU_Runner("./lib/")

# Load result calculated by the simplicial decomposition.
ndim, kdim, domain_div, hypercubes, vertices, simplices = runner.read_simplicial_decomp_binary("out/results.spx_bin")

# Select if you want to plot tetrahedra (True) or edges only (False) 
plot_tetrahedra = ((ndim - kdim) == 3)


# ==============================================================
# PyVista plotting functions
# ==============================================================
plotter = pv.Plotter()

if plot_tetrahedra == False:
    #----------------------------
    # Plotting the edges only 
    # ---------------------------
    # Preprocess/project vertices from n-dimension to 3D and create edges.
    # Also calculate color map. One value per vertex.
    segments = []
    vertex_cmap = []
    for i in range(simplices.shape[0]):
        for j in range(0, (ndim-kdim+1)-1):
            for k in range(j+1, (ndim-kdim+1)):
                vj = simplices[i,j]
                vk = simplices[i,k]
                pj = vertices[vj,:]
                pk = vertices[vk,:]
                segment = [pj[0], pj[1], pj[2], pk[0], pk[1], pk[2]]
                segments.append(segment)
                # vertex_cmap.append([pj[2], pk[2]])
                vertex_cmap.append([i, i])
                
    segments = np.array(segments)
    vertex_cmap = np.array(vertex_cmap)

    # Convert segments to lines format for PyVista.
    lines = []
    for i in range(segments.shape[0]):
        lines.extend([2, i*2, i*2+1])  # 2 points per line.

    points = segments.reshape(-1, 3)
    poly = pv.PolyData(points, lines=lines)
    bounds = poly.bounds

    # Setting color map.
    poly.point_data["vcmap"] = vertex_cmap.reshape((-1))

    # Add mesh and keep the actor so we can scale it.
    actor = plotter.add_mesh(
        poly,
        scalars="vcmap",            # uses point data.
        cmap="viridis",
        line_width=3,
        render_lines_as_tubes=False,
        scalar_bar_args={"title": "custom color map"}
    )
else:
    # ----------------------------
    # Plotting the tetrahedra
    # ----------------------------
    # Create the UnstructuredGrid object.
    # We pass the cell connectivity and cell types in a dictionary,
    # along with the points array.
    grid = pv.UnstructuredGrid({
            pv.CellType.TETRA: simplices[:,:4],
        },
        vertices[:,:3] # orthogonal projection to 3D
    )
    bounds = grid.bounds

    # Plot the resulting tetrahedral mesh with a clip plane.
    # Add mesh and keep the actor so we can scale it.
    actor = plotter.add_mesh_clip_plane(
        grid, 
        show_edges=True, # see the wireframe.
        crinkle=True, 
        interaction_event='always',
        normal=(0.,1.,0.)
    )

# ---- Matplotlib-like "auto" aspect (non-equal) via non-uniform scaling ----
xmin, xmax, ymin, ymax, zmin, zmax = bounds
rx = xmax - xmin
ry = ymax - ymin
rz = zmax - zmin

# Avoid division by zero if an extent is flat.
eps = 1e-12
rx = max(rx, eps)
ry = max(ry, eps)
rz = max(rz, eps)

# Scale smaller axes up so display extents match the largest axis.
L = max(rx, ry, rz)
sx, sy, sz = (L / rx, L / ry, L / rz)
actor.SetScale(sx, sy, sz)

# Keep the camera centered on the (unscaled) data bounds.
cx = 0.5 * (xmin + xmax)
cy = 0.5 * (ymin + ymax)
cz = 0.5 * (zmin + zmax)

plotter.camera.focal_point = (cx, cy, cz)
# A reasonable oblique/isometric-like position; tweak as you like.
plotter.camera.position = (cx + 1.5*L, cy + 1.5*L, cz + 1.5*L)
plotter.camera.up = (0, 0, 1)

# Reset after scaling so the scene fits the view.
plotter.renderer.reset_camera()

plotter.show()
