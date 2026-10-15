import os as _os
REPO = _os.environ.get('IMPLICIT_MANIFOLDS_ROOT') or next(_p for _p in (_os.path.abspath(_os.path.join(_os.path.dirname(__file__), *(['..'] * _i))) for _i in range(8)) if _os.path.exists(_os.path.join(_p, 'CITATION.cff')))
import sys, os
sys.path.insert(0, REPO + '/python/implicit_manifolds')
from render import *
from polfile import euler
R = REPO + '/book/runs/hypercubes/'
FG = REPO + '/book/figures/'
which = sys.argv[1:]
DONE = open(R+'timings.txt').read()
def ok(f): return os.path.basename(f) in DONE



def s3_figure(pol, out, cubes_sel):
    n, k, D, cells = read_pol(pol)
    fig = plt.figure(figsize=(10, 4.6))
    sel = [c for c in cells if tuple(c.grid) in cubes_sel]
    P = polygons(sel, dedupe=True)
    ax = new_ax(fig, 121, elev=18, azim=-50)
    pr = axes_proj(4, (0, 1, 2))
    draw_polys(ax, P, pr, color_values=[p[:, 3].mean() for p in P], cmap='coolwarm',
               alpha=0.55, lw=0.5, zoom=1.2)
    ax.set_title(f'cells of {len(set(tuple(c.grid) for c in sel))} grid cells ({len(sel)} polytopes)', fontsize=9)
    S = slice3(cells, 3, 0.3)
    ax = new_ax(fig, 122, elev=18, azim=-50)
    draw_polys(ax, S, pr, color_values=[p[:, 2].mean() for p in S], cmap='Blues', vmin=-2.5, vmax=1.2,
               lw=0.3, zoom=1.35)
    ax.set_title(f'slice $x_4=0.3$ ({len(S)} polygons)', fontsize=9)
    plt.subplots_adjust(0, 0, 1, 0.93, 0, 0)
    plt.savefig(FG + out, dpi=200)
    plt.close()


def stereo(P, R0):
    # stereographic projection of points of the 3-sphere of radius R0 in R^4 from (0,0,0,R0)
    return np.array([p[:3] * R0 / (R0 - p[3]) for p in P])


def surf4(pol, out, idx=(0, 1, 2), cvar=3, cmap='viridis', elev=25, azim=-60, title=None, ax=None, fig=None, zoom=1.3, lw=0.25, st=None):
    n, k, D, cells = read_pol(pol)
    P = polygons(cells)
    if st is not None:
        cv = [p[:, cvar].mean() for p in P]
        P = [stereo(p, st) for p in P]
        n = 3; idx = (0, 1, 2)
    own = ax is None
    if own:
        fig = plt.figure(figsize=(7, 5.2))
        ax = new_ax(fig, 111, elev=elev, azim=azim)
    draw_polys(ax, P, axes_proj(n, idx), color_values=(cv if st is not None else [p[:, cvar].mean() for p in P]), cmap=cmap, lw=lw, zoom=zoom)
    if title:
        ax.set_title(title, fontsize=9)
    if own:
        plt.subplots_adjust(0, 0, 1, 1)
        plt.savefig(FG + out, dpi=200)
        plt.close()
    return len(P)


if 's3' in which:
    n_, k_, D_, cc = read_pol(R + 'sph4_GMHC.pol')
    g0 = tuple(cc[len(cc)//3].grid)
    sel = {tuple(c.grid) for c in cc if max(abs(a-b) for a, b in zip(c.grid, g0)) <= 1 and c.grid[3] == g0[3]}
    sel = set(sorted(sel)[:6])
    C = REPO + '/book/runs/fch_book/BOOK_SPHERE4/pta.pol'
    if os.path.exists(C) and os.path.getsize(C) > 0:
        s3_figure(C, 's3_ms.png', sel)
    s3_figure(R + 'sph4_GMHC.pol', 's3_cmh.png', sel)

if 'torus' in which:
    for f, o in ((R + 'tor4_GCHC.pol', 'torus4_hypercube.png'), (REPO + '/book/runs/fch_book/BOOK_TORUS4/pta.pol', 'torus4_simplex.png')):
        if os.path.exists(f):
            m = surf4(f, o, cvar=2, elev=35, azim=-55, st=np.sqrt(1.25), cmap='viridis')
            print(o, m)

if 'complex' in which:
    CB = REPO + '/book/runs/fch_book/'
    for alg, out in (('fch', 'complex_gcch.png'), ('pta', 'complex_cs.png')):
        files = [CB + f'BOOK_{c}/{alg}.pol' for c in ('ZCOSW', 'Z2COSW', 'Z2W5')]
        if not all(os.path.exists(f) and os.path.getsize(f) > 0 for f in files):
            continue
        fig = plt.figure(figsize=(11, 4))
        for i, f in enumerate(files):
            ax = new_ax(fig, 131 + i, elev=22, azim=-60)
            m = surf4(f, None, idx=(0, 1, 2), cvar=3, ax=ax, fig=fig, zoom=1.25, lw=0.15)
            print(out, i, m)
        plt.subplots_adjust(0, 0, 1, 1, 0, 0)
        plt.savefig(FG + out, dpi=200)
        plt.close()

if 'klein' in which:
    n, k, D, cells = read_pol(REPO + '/book/runs/fch/KLEIN_5TO3/fch.pol')
    P = polygons(cells)
    fig = plt.figure(figsize=(7, 5))
    ax = new_ax(fig, 111, elev=28, azim=-35)
    draw_polys(ax, P, axes_proj(5, (0, 1, 2)), color_values=[p[:, 4].mean() for p in P], cmap='twilight',
               alpha=0.92, lw=0.12, zoom=1.35)
    plt.subplots_adjust(0, 0, 1, 1)
    plt.savefig(FG + 'klein.png', dpi=200)
    plt.close()
    print('klein', len(P))

if 'springs' in which:
    B = REPO + '/python/fch_lib/out/'
    fig = plt.figure(figsize=(11, 4.4))
    n, k, D, cells = read_pol(B + 'spring1d/results.pol')
    S = segments(cells)
    ax = new_ax(fig, 131, elev=15, azim=-60)
    draw_segments(ax, S, axes_proj(3), color='C3', lw=1.6)
    n, k, D, cells = read_pol(B + 'spring2d/results.pol')
    P = polygons(cells)
    ax = new_ax(fig, 132, elev=15, azim=-60)
    draw_polys(ax, P, axes_proj(3), color_values=[p[:, 2].mean() for p in P], cmap='viridis', lw=0.15, zoom=1.2)
    n, k, D, cells = read_pol(B + 'spring3d/results.pol')
    S3 = slice3(cells, 3, 0.9)
    ax = new_ax(fig, 133, elev=15, azim=-60)
    draw_polys(ax, S3, axes_proj(4, (0, 1, 2)), color_values=[p[:, 2].mean() for p in S3], cmap='plasma', lw=0.12, zoom=1.2)
    plt.subplots_adjust(0, 0, 1, 1, 0, 0)
    plt.savefig(FG + 'springs.png', dpi=200)
    plt.close()

if 'solid' in which:
    n, k, D, cells = read_pol(REPO + '/python/fch_lib/out/solidspheres/results.pol')
    fig = plt.figure(figsize=(11, 3.8))
    for i, c in enumerate((0.85, 0.95, 0.999)):
        S = slice3(cells, 3, c)
        ax = new_ax(fig, 131 + i, elev=25, azim=-60)
        draw_polys(ax, S, axes_proj(4, (0, 1, 2)), color_values=[p[:, 2].mean() for p in S], cmap='Blues', vmin=-2, vmax=1.2, lw=0.2, zoom=1.3)
        ax.set_title(f'$x_4={c:.2f}$', fontsize=9)
        ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_zlim(-2, 2)
    plt.subplots_adjust(0, 0, 1, 0.92, 0, 0)
    plt.savefig(FG + 'solidspheres.png', dpi=200)
    plt.close()


def pair(files, out, titles, elev=25, azim=-60, size=(10, 4.4), cmap='Blues', cvar=2):
    fig = plt.figure(figsize=size)
    for i, f in enumerate(files):
        n, k, D, cells = read_pol(R + f)
        P = polygons(cells)
        ax = new_ax(fig, 100 + 10 * len(files) + 1 + i, elev=elev, azim=azim)
        allz = np.concatenate([p[:, cvar] for p in P])
        draw_polys(ax, P, axes_proj(3), color_values=[p[:, cvar].mean() for p in P], cmap=cmap,
                   vmin=allz.min() - 0.8 * np.ptp(allz), vmax=allz.max(), lw=0.15, zoom=1.3)
        ax.set_title(titles[i] + f' ({len(P)} cells)', fontsize=9)
        print(out, f, len(P), euler(cells)[1])
    plt.subplots_adjust(0, 0, 1, 0.93, 0, 0)
    plt.savefig(FG + out, dpi=200)
    plt.close()


if 'model' in which:
    if ok('two_p8_GMHC.pol'):
        pair(['two_min_MS.pol', 'two_min_GMHC.pol', 'two_p8_GMHC.pol'], 'two_spheres.png',
             ['union, simplices', 'union, hypercubes', '$p=8$, hypercubes'], size=(11, 3.8))
    if ok('toro8_GMHC.pol'):
        pair(['toro8_MS.pol', 'toro8_GMHC.pol'], 'eight_tori.png', ['simplices', 'hypercubes'], elev=40)
    if ok('knot_GMHC.pol'):
        pair(['knot_MS.pol', 'knot_GMHC.pol'], 'knot.png', ['simplices', 'hypercubes'], elev=35)
    if ok('intesf_GMHC.pol'):
        pair(['intesf_MS.pol', 'intesf_GMHC.pol'], 'intersection.png', ['simplices', 'hypercubes'], elev=20)
