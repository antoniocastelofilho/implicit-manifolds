import numpy as np
from fch_wrapper import FCH_CPU_Runner
import pyvista as pv

runner = FCH_CPU_Runner("./lib/")

# Load result calculated by fch.
edges = runner.read_fch_output_binary("out/results.bin")

# Preprocess/project edges from n-dimension to 3D.
ndim = int(edges.shape[1] / 2)
segments = np.concatenate((edges[:,0:3], edges[:,ndim+0:ndim+3]), axis=1)

# Calculate color map. One value per vertex.
vertex_cmap = np.concatenate((edges[:,2].reshape(-1,1), edges[:,ndim+2].reshape(-1,1)), axis=1)

# ==============================================================
# PyVista plotting functions
# ==============================================================
# Convert segments to lines format for PyVista.
lines = []
for i in range(segments.shape[0]):
    lines.extend([2, i*2, i*2+1])  # 2 points per line.

points = segments.reshape(-1, 3)
poly = pv.PolyData(points, lines=lines)

# Setting color map.
poly.point_data["vcmap"] = vertex_cmap.reshape((-1))

plotter = pv.Plotter()

# Add mesh and keep the actor so we can scale it.
actor = plotter.add_mesh(
    poly,
    scalars="vcmap",            # uses point data.
    cmap="viridis",
    line_width=3,
    render_lines_as_tubes=False,
    scalar_bar_args={"title": "custom color map"}
)

# ---- Matplotlib-like "auto" aspect (non-equal) via non-uniform scaling ----
xmin, xmax, ymin, ymax, zmin, zmax = poly.bounds
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
