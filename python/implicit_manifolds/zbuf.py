"""Small z-buffer renderer (orthographic) for polygon soups, with flat two-sided shading,
depth-tested cell outlines and supersampling.  Avoids the sorting artefacts of matplotlib 3D."""
import numpy as np
from numba import njit

def view_matrix(elev, azim):
    e, a = np.radians(elev), np.radians(azim)
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(np.pi/2 - e), -np.sin(np.pi/2 - e)], [0, np.sin(np.pi/2 - e), np.cos(np.pi/2 - e)]])
    # camera looks along -z after rotation: x right, y up
    return Rx @ Rz

@njit(cache=True)
def _raster(T, C, W, H, zb, img):
    for t in range(T.shape[0]):
        x0, y0, z0 = T[t, 0]; x1, y1, z1 = T[t, 1]; x2, y2, z2 = T[t, 2]
        minx = max(int(np.floor(min(x0, x1, x2))), 0); maxx = min(int(np.ceil(max(x0, x1, x2))), W - 1)
        miny = max(int(np.floor(min(y0, y1, y2))), 0); maxy = min(int(np.ceil(max(y0, y1, y2))), H - 1)
        den = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        if abs(den) < 1e-12:
            continue
        for py in range(miny, maxy + 1):
            for px in range(minx, maxx + 1):
                fx = px + 0.5; fy = py + 0.5
                l0 = ((y1 - y2) * (fx - x2) + (x2 - x1) * (fy - y2)) / den
                l1 = ((y2 - y0) * (fx - x2) + (x0 - x2) * (fy - y2)) / den
                l2 = 1.0 - l0 - l1
                if l0 < 0 or l1 < 0 or l2 < 0:
                    continue
                z = l0 * z0 + l1 * z1 + l2 * z2
                if z > zb[py, px]:
                    zb[py, px] = z
                    img[py, px, 0] = C[t, 0]; img[py, px, 1] = C[t, 1]; img[py, px, 2] = C[t, 2]

@njit(cache=True)
def _lines(L, W, H, zb, img, col, width, bias):
    for s in range(L.shape[0]):
        x0, y0, z0 = L[s, 0]; x1, y1, z1 = L[s, 1]
        n = int(max(abs(x1 - x0), abs(y1 - y0)) * 2) + 2
        for i in range(n + 1):
            t = i / n
            x = x0 + t * (x1 - x0); y = y0 + t * (y1 - y0); z = z0 + t * (z1 - z0)
            for dx in range(-width, width + 1):
                for dy in range(-width, width + 1):
                    px = int(x) + dx; py = int(y) + dy
                    if 0 <= px < W and 0 <= py < H and z + bias >= zb[py, px]:
                        img[py, px, 0] = col[0]; img[py, px, 1] = col[1]; img[py, px, 2] = col[2]

def render(polys, colors, elev=25, azim=-60, size=(1200, 900), ss=2, edges=True, edge_color=(0.1, 0.1, 0.1),
           edge_alpha=0.35, light=(0.3, -0.4, 0.85), ambient=0.38, margin=0.04, bg=(1, 1, 1), R=None, lw=0):
    if R is None:
        R = view_matrix(elev, azim)
    W, H = size[0] * ss, size[1] * ss
    P3 = [np.asarray(p, float) @ R.T for p in polys]
    allp = np.vstack(P3); lo = allp.min(0); hi = allp.max(0)
    scale = (1 - 2 * margin) * min(W / (hi[0] - lo[0] + 1e-12), H / (hi[1] - lo[1] + 1e-12))
    cx = (lo[0] + hi[0]) / 2; cy = (lo[1] + hi[1]) / 2
    def scr(Q):
        return np.c_[(Q[:, 0] - cx) * scale + W / 2, H / 2 - (Q[:, 1] - cy) * scale, Q[:, 2] * scale]
    Lv = np.asarray(light, float); Lv /= np.linalg.norm(Lv)
    tris = []; cols = []; segs = []
    for Q, c in zip(P3, colors):
        S = scr(Q); g = Q.mean(0)
        n = np.zeros(3)
        for i in range(len(Q)):
            n += np.cross(Q[i] - g, Q[(i + 1) % len(Q)] - g)
        nn = np.linalg.norm(n)
        sh = ambient + (1 - ambient) * (abs(n @ Lv) / nn if nn > 0 else 0.5)
        col = np.clip(np.asarray(c[:3]) * sh, 0, 1)
        gs = np.array([(g[0] - cx) * scale + W / 2, H / 2 - (g[1] - cy) * scale, g[2] * scale])
        for i in range(len(S)):
            tris.append([gs, S[i], S[(i + 1) % len(S)]]); cols.append(col)
            if edges: segs.append([S[i], S[(i + 1) % len(S)]])
    T = np.array(tris, float); C = np.array(cols, float)
    zb = np.full((H, W), -1e300); img = np.ones((H, W, 3)) * np.asarray(bg)
    _raster(T, C, W, H, zb, img)
    if edges and segs:
        base = img.copy()
        _lines(np.array(segs, float), W, H, zb, img, np.asarray(edge_color, float), lw, 1.5 * ss)
        img = base * (1 - edge_alpha) + img * edge_alpha
    # downsample
    img = img.reshape(size[1], ss, size[0], ss, 3).mean(axis=(1, 3))
    return img

def crop(img, pad=10, bg=1.0):
    m = (np.abs(img - bg).sum(2) > 0.02)
    ys, xs = np.where(m)
    if len(xs) == 0: return img
    return img[max(ys.min() - pad, 0):ys.max() + pad, max(xs.min() - pad, 0):xs.max() + pad]
