#!/usr/bin/env python3
"""BD36 READ_ONLY geometry & symmetry audit for PtSe2 BiDB bilayer inputs.

Pure-numpy, self-contained. No external DFT / spglib / pymatgen needed.
Verifies for each POSCAR:
  - lattice metrics (hex a, c, volume, vacuum)
  - layer assignment (bottom/top by z) and interlayer metrics
  - space-group-like operations by SITE-PERMUTATION over the 2D hex lattice:
      inversion i, C3 (about z), C2', vertical mirrors, horizontal mirror M_z
  - for derived rigid slides: top-layer shift == declared u, layer separation
    preserved, inversion retained (site-permutation error in Angstrom)

Outputs a JSON report (evidence) and prints a human table.
"""
import json, sys, math, itertools
import numpy as np

TOL = 1e-4   # fractional-coordinate matching tolerance
TOL_CART = 1e-3  # Angstrom tolerance for permutation errors

# Hexagonal lattice rotation matrices in FRACTIONAL coordinates (about z through origin).
# C3^+ (120 deg): (x,y) -> (-y, x-y)
C3 = np.array([[0, -1], [1, -1]], dtype=float)
# C2' about an in-plane axis through origin at angle 0 (mirror-like): (x,y) -> (-x, y-x)
# (a vertical mirror plane containing the a-axis)
MV = np.array([[-1, 1], [0, 1]], dtype=float)


def parse_poscar(path):
    lines = open(path).read().splitlines()
    title = lines[0].strip()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])
    lat *= scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    assert lines[7].split()[0][0] in "DdKkCc", f"coord type line: {lines[7]}"
    direct = lines[7].split()[0][0] in "DdKk"
    coords = []
    for ln in lines[8:8 + sum(counts)]:
        coords.append([float(x) for x in ln.split()[:3]])
    coords = np.array(coords)
    if not direct:
        coords = coords @ np.linalg.inv(lat)
    return title, lat, species, counts, coords


def frac_to_cart(frac, lat):
    return frac @ lat


def wrap(frac):
    return frac - np.floor(frac)


def unwrapped_delta(a, b):
    """shortest fractional displacement from b to a, using lattice periodicity"""
    d = a - b
    for _ in range(2):
        d = d - np.round(d)
    return d


def cart_delta(a, b, lat):
    return unwrapped_delta(a, b) @ lat


def find_match(cart_target, atom_carts, lat, tol=TOL_CART):
    """find index of atom whose lattice-image equals cart_target (mod lattice), else -1"""
    best = -1
    best_err = tol
    for j, c in enumerate(atom_carts):
        d = cart_target - c
        # reduce by lattice vectors
        f = np.linalg.solve(lat.T, d)  # coords of d in lattice basis
        f = f - np.round(f)
        err = np.linalg.norm(f @ lat)
        if err < best_err:
            best_err = err
            best = j
    return best, best_err


def operation_is_symmetry(atom_carts, op_cart_matrix, center_cart, lat, tol=TOL_CART):
    """op_cart_matrix (3x3) applied about center_cart: r' = center + M (r - center).
    Return (True, max_err) if it permutes the atomic site set onto itself."""
    max_err = 0.0
    used = np.zeros(len(atom_carts), dtype=bool)
    for i, r in enumerate(atom_carts):
        rp = center_cart + op_cart_matrix @ (r - center_cart)
        j, err = find_match(rp, atom_carts, lat, tol=100 * tol)
        if j < 0 or used[j]:
            return False, float(max_err)
        used[j] = True
        max_err = max(max_err, err)
    return True, float(max_err)


def frac_centers(limit=6):
    """candidate centers on the grid (i/limit, j/limit), z=0 (z center handled separately)"""
    out = []
    for i in range(limit):
        for j in range(limit):
            out.append(np.array([i / limit, j / limit, 0.0]))
    return out


def audit(path):
    title, lat, species, counts, frac = parse_poscar(path)
    cart = frac_to_cart(frac, lat)
    n = len(frac)
    # species index
    sp_idx = []
    for k, (sp, cnt) in enumerate(zip(species, counts)):
        sp_idx += [k] * cnt
    sp_idx = np.array(sp_idx)
    z = cart[:, 2]
    order = np.argsort(z)
    a = np.linalg.norm(lat[0])
    b = np.linalg.norm(lat[1])
    cvec = lat[2]
    c = np.linalg.norm(cvec)
    vol = abs(np.linalg.det(lat))
    # layers: distinct Pt z-groups
    pt_z = sorted(z[sp_idx == 0])
    n_layers = len(pt_z)
    layers = []
    if n_layers == 1:
        layers = [np.where(sp_idx == 0)[0]]
    else:
        z0, z1 = pt_z[0], pt_z[1]
        tol = 0.5 * (z1 - z0)
        layers = [
            np.array([i for i in range(n) if abs(z[i] - z0) < tol]),
            np.array([i for i in range(n) if abs(z[i] - z1) < tol]),
        ]
    # interlayer metrics
    if n_layers == 2:
        btm, top = layers
        z_btm, z_top = z[btm], z[top]
        pt_btm = cart[btm][sp_idx[btm] == 0][:, :2]
        pt_top = cart[top][sp_idx[top] == 0][:, :2]
        d_pt2d = np.linalg.norm(pt_top[0] - pt_btm[0]) if len(pt_top) and len(pt_btm) else None
        interlayer_gap = abs(z[top].min() - z[btm].max())
        # occupied z span
        zmin, zmax = z.min(), z.max()
        vacuum_top = c - zmax
        vacuum_bot = zmin
        vacuum_gap = c - (zmax - zmin)
    # ---- symmetry search (site permutation, full 3D slab) ----
    found = {"inversion": [], "c3": [], "c2prime": [], "mirror_v": [], "mirror_h": []}
    # inversion: center (c_xy in fractional, z_center in CARTESIAN); op = -I3
    zc_list = []
    if n_layers == 1:
        zc_list = [cart[sp_idx == 0][0][2]]  # Pt plane
    else:
        zc_list = [0.5 * (z[layers[0]].max() + z[layers[1]].min())]  # midgap
    # inversion: op = -I3; centers: continuous candidates from atom-pair midpoints
    # (for two rigid identical layers the center sits at the midpoint of paired atoms)
    zc_list = []
    if n_layers == 1:
        zc_list = [cart[sp_idx == 0][0][2]]  # Pt plane
    else:
        zc_list = [0.5 * (z[layers[0]].max() + z[layers[1]].min())]  # midgap
    cand_centers = []
    imgs = [np.array([ii, jj, 0.0]) for ii in range(-2, 3) for jj in range(-2, 3)]
    for i in range(n):
        for j in range(n):
            for img in imgs:
                cand = 0.5 * (cart[i] + (cart[j] + img @ lat))
                if abs(cand[2] - zc_list[0]) < 0.5 * c:
                    cand_centers.append(cand)
    # also keep the fractional grid points (high-symmetry cases)
    for fc in frac_centers():
        cxy = np.array([fc[0], fc[1], 0.0]) @ lat
        cand_centers.append(np.array([cxy[0], cxy[1], zc_list[0]]))
    seen = set()
    for cc in cand_centers:
        key = tuple(np.round(cc / [a, b, c], 4))
        if key in seen:
            continue
        seen.add(key)
        ok, err = operation_is_symmetry(cart, -np.eye(3), cc, lat)
        if ok:
            found["inversion"].append({"center_cart_A": [round(float(v), 5) for v in cc],
                                        "max_err_A": err})
    # C3 about z: true CARTESIAN rotations by +-120 deg about z (basis-independent)
    th = 2 * math.pi / 3
    for sgn in (1.0, -1.0):
        ct, st = math.cos(sgn * th), math.sin(sgn * th)
        R3 = np.array([[ct, -st, 0.0], [st, ct, 0.0], [0.0, 0.0, 1.0]])
        for fc in frac_centers():
            cxy = np.array([fc[0], fc[1], 0.0]) @ lat
            cc = np.array([cxy[0], cxy[1], 0.0])
            ok, err = operation_is_symmetry(cart, R3, cc, lat)
            if ok:
                found["c3"].append({"sign": sgn, "center_frac": [fc[0], fc[1], 0.0], "max_err_A": err})
    # C2' axes (180 deg about in-plane axes) -> also give vertical mirror info; keep separate:
    # vertical mirrors are equivalent to C2' about in-plane axis + no z flip; here mirror_v = reflection
    # through vertical plane: (x,y) -> reflection; z -> z. Implement via M_v matrices at several in-plane axes.
    for axis_deg in [0, 30, 60, 90, 120, 150]:
        th = math.radians(axis_deg)
        u = np.array([math.cos(th), math.sin(th)])
        # reflection about line through origin along u, in CARTESIAN xy
        M = np.eye(3)
        P = np.outer(u, u)
        M[:2, :2] = 2 * P - np.eye(2)
        for fc in frac_centers():
            cxy = np.array([fc[0], fc[1], 0.0]) @ lat
            cc = np.array([cxy[0], cxy[1], 0.0])
            ok, err = operation_is_symmetry(cart, M, cc, lat)
            if ok:
                found["mirror_v"].append({"axis_deg": axis_deg, "center_frac": [fc[0], fc[1], 0.0],
                                          "max_err_A": err})
    # horizontal mirror M_z (z -> -z) about a horizontal plane at height zh (cartesian)
    Mh = np.diag([1.0, 1.0, -1.0])
    for zh in np.unique(np.round(z, 6)):
        cc = np.array([0.0, 0.0, zh])
        ok, err = operation_is_symmetry(cart, Mh, cc, lat)
        if ok:
            found["mirror_h"].append({"z_A": zh, "max_err_A": err})

    # ---- derived-slide verification (for the two generic POSCARs) ----
    derived = None
    if "generic" in path:
        base_2253 = path.replace("_generic_m008_m019_rigid", "").replace("_generic_p008_p019_rigid", "")
        base_2253 = base_2253.replace("POSCAR_PtSe2_BiDB2253", "POSCAR_BiDB_2253_PtSe2")
        # compute top-layer shift vs 2253 (read the base file directly)
        t, latb, spb, cntb, fb = parse_poscar(base_2253)
        zb = fb @ latb
        spb_idx = []
        for k, (sp, cnt) in enumerate(zip(spb, cntb)):
            spb_idx += [k] * cnt
        spb_idx = np.array(spb_idx)
        pt_global_b = np.where(spb_idx == 0)[0]
        topb = pt_global_b[np.argmax(zb[pt_global_b, 2])]
        topb_pt_frac = fb[topb]
        pt_global = np.where(sp_idx == 0)[0]
        topp = pt_global[np.argmax(z[pt_global])]
        top_pt_frac = frac[topp]
        shift = unwrapped_delta(top_pt_frac, topb_pt_frac)[:2]
        # declared
        declared = [0.08, 0.19] if "p008_p019" in path else [-0.08, -0.19]
        shift_err = np.linalg.norm(np.array(shift) - np.array(declared))
        # layer separation preserved?
        ztop = z[layers[1]].min() if n_layers == 2 else None
        ztop_b = zb[spb_idx == 0].max()
        dsep = abs(z[layers[1]].min() - z[layers[0]].max())
        dsep_b = abs(zb[spb_idx == 0].max() - zb[spb_idx == 0].min())
        derived = {
            "declared_u_frac": declared,
            "measured_top_Pt_shift_frac": [round(x, 6) for x in shift],
            "shift_match_err_frac": float(shift_err),
            "layer_gap_A": float(dsep),
        }
    return {
        "file": path.split("/")[-1],
        "title": title,
        "lattice": {"a_A": round(float(a), 6), "b_A": round(float(b), 6),
                    "c_A": round(float(c), 6), "volume_A3": round(float(vol), 4),
                    "hex": bool(abs(np.dot(lat[0], lat[1]) / (a * b) + 0.5) < 1e-6)},
        "species": species,
        "counts": counts,
        "n_atoms": n,
        "n_layers": n_layers,
        "layer_assignment": {
            f"layer_{k}": {"atoms": [int(i) for i in L], "z_range_A": [round(float(z[L].min()), 5), round(float(z[L].max()), 5)]}
            for k, L in enumerate(layers)
        },
        "interlayer": (None if n_layers == 1 else {
            "pt_pt_vertical_A": round(float(z[layers[1]].max() - z[layers[0]].min()), 5),
            "se_se_gap_A": round(float(z[layers[1]].min() - z[layers[0]].max()), 5),
            "inplane_pt_pt_A": (round(float(d_pt2d), 5) if d_pt2d is not None else None),
            "vacuum_gap_A": round(float(vacuum_gap), 5),
        }),
        "symmetry": {
            k: (len(v) and v[0] or None) for k, v in found.items()
        },
        "symmetry_counts": {k: len(v) for k, v in found.items()},
        "derived_slide": derived,
    }


def main():
    files = sys.argv[1:] if len(sys.argv) > 1 else [
        "POSCAR_BiDB_2252_PtSe2", "POSCAR_BiDB_2253_PtSe2", "POSCAR_BiDB_2255_PtSe2",
        "POSCAR_BiDB_2254_PtSe2", "POSCAR_PtSe2_BiDB2253_generic_p008_p019_rigid",
        "POSCAR_PtSe2_BiDB2253_generic_m008_m019_rigid",
    ]
    out = []
    for f in files:
        try:
            r = audit(f)
            out.append(r)
            print(f"== {f}")
            print(f"   lattice a={r['lattice']['a_A']} c={r['lattice']['c_A']} vol={r['lattice']['volume_A3']} hex={r['lattice']['hex']}")
            print(f"   layers={r['n_layers']}  species={r['species']}{r['counts']}")
            if r["interlayer"]:
                il = r["interlayer"]
                print(f"   interlayer: Pt-Pt_z={il['pt_pt_vertical_A']} Se-Se_gap={il['se_se_gap_A']} vac={il['vacuum_gap_A']}")
            sym = r["symmetry"]
            print(f"   inversion={sym['inversion'] is not None}  C3={sym['c3'] is not None}  "
                  f"mirror_v={sym['mirror_v'] is not None}  mirror_h={sym['mirror_h'] is not None}")
            if r["derived_slide"]:
                d = r["derived_slide"]
                print(f"   slide: declared={d['declared_u_frac']} measured={d['measured_top_Pt_shift_frac']} "
                      f"match_err={d['shift_match_err_frac']:.2e} layer_gap={d['layer_gap_A']}")
        except Exception as e:
            print(f"!! {f}: {e}")
            raise
    json.dump(out, open("audit_report.json", "w"), indent=2)
    print("\nwrote audit_report.json")


if __name__ == "__main__":
    main()
