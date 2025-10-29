#!/usr/bin/env python3
"""
Generate Lebedev C++ headers and an aggregator source file, including face indices.

- Outputs:
    <OUT_DIR>/cpp/include/sasktran2/math/unitsphere/lebedev/sphere_<N>.h
    <OUT_DIR>/cpp/lib/unitsphere/lebedev_autogen.cpp

- Dependencies:
    pip install numpy scipy quadpy
"""

import argparse
import math
from pathlib import Path
from typing import List, Tuple

import numpy as np
from scipy.integrate import lebedev_rule
from scipy.spatial import ConvexHull


# Default sizes: matches your set
DEFAULT_NS = [
    6,
    14,
    26,
    38,
    50,
    74,
    86,
    110,
    146,
    170,
    194,
    230,
    266,
    302,
    350,
    434,
    590,
    770,
    974,
    1202,
    1454,
    1730,
    2030,
    2354,
    2702,
    3074,
    3470,
    3890,
]

# SciPy's lebedev_rule expects an odd "degree" and returns N implied by that degree.
SP_MS = [
    6,
    14,
    26,
    38,
    50,
    74,
    86,
    110,
    146,
    170,
    194,
    230,
    266,
    302,
    350,
    434,
    590,
    770,
    974,
    1202,
    1454,
    1730,
    2030,
    2354,
    2702,
    3074,
    3470,
    3890,
    4334,
    4802,
    5294,
    5810,
]
SP_NS = [
    3,
    5,
    7,
    9,
    11,
    13,
    15,
    17,
    19,
    21,
    23,
    25,
    27,
    29,
    31,
    35,
    41,
    47,
    53,
    59,
    65,
    71,
    77,
    83,
    89,
    95,
    101,
    107,
    113,
    119,
    125,
    131,
]
SP_MAP = dict(zip(SP_MS, SP_NS))

# C++ namespace & include layout (adjust if yours differs)
NS = "sasktran2::math::unitsphere::lebedev"
INCLUDE_PREFIX = "sasktran2/math/unitsphere/lebedev"
INTERNAL_COMMON = "<sasktran2/internal_common.h>"
LEBEDEV_PUBLIC_H = "<sasktran2/math/unitsphere/lebedev.h>"


def c_double(x: float) -> str:
    s = f"{float(x):.17g}"
    if ("e" not in s) and ("E" not in s) and ("." not in s):
        s += ".0"
    return s.replace("E", "e")


def format_flat_array(vals: List[str], indent: int = 53, per_line: int = 4) -> str:
    pad = " " * indent
    out_lines = []
    for i in range(0, len(vals), per_line):
        chunk = ", ".join(vals[i : i + per_line])
        out_lines.append(pad + chunk + ",")
    if out_lines:
        out_lines[-1] = out_lines[-1].rstrip(",")
    return "\n".join(out_lines)


def format_xyzw_block(xyz: np.ndarray, w: np.ndarray) -> str:
    n = xyz.shape[1]
    vals: List[str] = []
    for i in range(n):
        vals.append(c_double(xyz[0, i]))
        vals.append(c_double(xyz[1, i]))
        vals.append(c_double(xyz[2, i]))
        vals.append(c_double(w[i]))
    return format_flat_array(vals, indent=53, per_line=4)


def format_faces_block(faces: np.ndarray) -> str:
    # faces shape: (F, 3), int
    vals = [str(int(v)) for v in faces.reshape(-1)]
    # use more per-line items for compactness (e.g., 12 = 4 triangles/line)
    return format_flat_array(vals, indent=53, per_line=12)


def orient_faces_ccw_outward(xyz: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """
    Ensure each row [i,j,k] is CCW with outward normal, i.e., (a x b) · c > 0.
    xyz: (3, N), faces: (F,3)
    """
    A = xyz[:, faces[:, 0]].T  # (F,3)
    B = xyz[:, faces[:, 1]].T
    C = xyz[:, faces[:, 2]].T
    cross = np.cross(A, B)  # (F,3)
    triple = np.einsum("ij,ij->i", cross, C)  # (F,)
    flip = triple < 0.0
    faces_flipped = faces.copy()
    faces_flipped[flip, 1], faces_flipped[flip, 2] = faces[flip, 2], faces[flip, 1]
    return faces_flipped


def compute_faces_from_convex_hull(xyz: np.ndarray) -> np.ndarray:
    """
    xyz: (3, N). Returns faces (F, 3) with CCW/outward orientation.
    """
    pts = xyz.T.copy()  # (N,3)
    hull = ConvexHull(pts, qhull_options="Qt")  # 'Qt' for triangulated facets
    faces = hull.simplices.copy()  # (F,3) int indices
    faces = orient_faces_ccw_outward(xyz, faces)
    return faces


def write_sphere_header(
    out_dir: Path, N: int, xyz: np.ndarray, w: np.ndarray, faces: np.ndarray
) -> Path:
    out_path = out_dir / f"cpp/include/sasktran2/math/unitsphere/lebedev/sphere_{N}.h"
    xyzw_body = format_xyzw_block(xyz, w)
    faces_body = format_faces_block(faces)
    F = faces.shape[0]

    text = f"""#pragma once

#include {INTERNAL_COMMON}

namespace {NS} {{
    // Packed as x0, y0, z0, w0, x1, y1, z1, w1, ...
    const static double g_lebedev_xyzw_{N}[4 * {N}] = {{
{xyzw_body}
    }};

    // Triangulated faces (CCW, outward), 3 * g_lebedev_num_faces_{N} entries
    const static int g_lebedev_num_faces_{N} = {F};
    const static int g_lebedev_faces_{N}[3 * g_lebedev_num_faces_{N}] = {{
{faces_body}
    }};
}} // namespace {NS}
"""
    out_path.write_text(text)
    return out_path


def write_aggregator(out_dir: Path, Ns: List[int]) -> Path:
    includes = "\n".join(f"#include <{INCLUDE_PREFIX}/sphere_{N}.h>" for N in Ns)

    chain_parts = []
    for i, N in enumerate(Ns):
        head = "if" if i == 0 else "else if"
        chain_parts.append(
            f"""        {head} (npoints == {N}) {{
            result = Eigen::Map<const Eigen::MatrixXd>(g_lebedev_xyzw_{N}, 4, {N});
            faces = Eigen::Map<const Eigen::MatrixXi>(g_lebedev_faces_{N}, 3, g_lebedev_num_faces_{N});
        }}"""
        )
    chain = "\n".join(chain_parts)

    text = f"""#include {INTERNAL_COMMON}
#include {LEBEDEV_PUBLIC_H}

{includes}

namespace {NS} {{
    void get_lebedev_data(int npoints, Eigen::MatrixXd& result, Eigen::MatrixXi& faces) {{
{chain}
        else {{
            spdlog::error(
                "Requested number of Lebedev quadrature points does not exist");
            throw std::runtime_error(
                "Requested number of Lebedev quadrature points does not exist");
        }}
    }}
}} // namespace {NS}
"""
    out_path = out_dir / "cpp/lib/unitsphere/lebedev_autogen.cpp"
    out_path.write_text(text)
    return out_path


def get_lebedev_xyz_w(N: int) -> Tuple[np.ndarray, np.ndarray]:
    # SciPy returns (xyz(3,N), w(N,)) for degree SP_MAP[N]
    xyz, w = lebedev_rule(SP_MAP[N])
    # normalize points (safety)
    norms = np.linalg.norm(xyz, axis=0)
    xyz[:, norms > 0] /= norms[norms > 0]
    # Normalize weights to integrate a function over the unit sphere:
    # SciPy's lebedev_rule integrates over the sphere surface area; divide by 4π
    return xyz, w / (4.0 * math.pi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-o",
        "--out-dir",
        type=str,
        required=True,
        help="Directory to write generated headers/sources into",
    )
    ap.add_argument(
        "--sizes",
        type=int,
        nargs="*",
        default=DEFAULT_NS,
        help="Lebedev sizes to generate",
    )
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    (out_dir / "cpp/include/sasktran2/math/unitsphere/lebedev").mkdir(
        parents=True, exist_ok=True
    )
    (out_dir / "cpp/lib/unitsphere").mkdir(parents=True, exist_ok=True)

    generated = []
    for N in args.sizes:
        xyz, w = get_lebedev_xyz_w(N)
        faces = compute_faces_from_convex_hull(xyz)  # (F,3), CCW outward
        path = write_sphere_header(out_dir, N, xyz, w, faces)
        generated.append((N, path))

    Ns_sorted = sorted(N for N, _ in generated)
    agg_path = write_aggregator(out_dir, Ns_sorted)

    print(f"Generated {len(generated)} headers and aggregator:")
    for N, p in sorted(generated):
        print(f"  - {p}")
    print(f"  - {agg_path}")


if __name__ == "__main__":
    main()
