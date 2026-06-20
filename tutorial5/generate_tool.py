#!/usr/bin/env python3
# Python Script
# Licensed under CC BY 4.0
# A Set of Materials-Science Tutorials for the LAMMPS Simulation Package
# Tutorial 6: write an STL triangle mesh for a "winking smiley" stamping tool:
# a round outer ridge framing one open eye, one winking eye, and a smile.  The
# features are pillars that point down (-z); LAMMPS create_atoms reads the mesh
# and fills its surface with (carbon) atoms.  Pressed into a metal surface, the
# tool leaves an emoji imprint.  ANY watertight STL works here -- e.g. the many
# models in online 3D-printing libraries -- so feel free to stamp your own.
#
#   python3 generate_tool.py            # writes tool.stl

import math

H = 12.0       # pillar height (z = 0 stamps the surface, z = H is the top)
nseg = 48      # segments around circular features


def facet(f, a, b, c):
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    m = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
    f.write(f"  facet normal {nx/m:.5e} {ny/m:.5e} {nz/m:.5e}\n")
    f.write("    outer loop\n")
    for v in (a, b, c):
        f.write(f"      vertex {v[0]:.5e} {v[1]:.5e} {v[2]:.5e}\n")
    f.write("    endloop\n  endfacet\n")


def cylinder(f, cx, cy, r, z0=0.0, z1=H, n=nseg):
    pts = [(cx + r * math.cos(2 * math.pi * i / n),
            cy + r * math.sin(2 * math.pi * i / n)) for i in range(n)]
    for i in range(n):
        j = (i + 1) % n
        a, b = pts[i], pts[j]
        facet(f, (a[0], a[1], z0), (b[0], b[1], z0), (b[0], b[1], z1))
        facet(f, (a[0], a[1], z0), (b[0], b[1], z1), (a[0], a[1], z1))
        facet(f, (cx, cy, z0), (b[0], b[1], z0), (a[0], a[1], z0))
        facet(f, (cx, cy, z1), (a[0], a[1], z1), (b[0], b[1], z1))


def ring(f, cx, cy, ro, ri, z0=0.0, z1=H, n=nseg):
    out = [(cx + ro * math.cos(2 * math.pi * i / n),
            cy + ro * math.sin(2 * math.pi * i / n)) for i in range(n)]
    inn = [(cx + ri * math.cos(2 * math.pi * i / n),
            cy + ri * math.sin(2 * math.pi * i / n)) for i in range(n)]
    for i in range(n):
        j = (i + 1) % n
        for (p, q, zlo, zhi) in [(out[i], out[j], z0, z1), (inn[j], inn[i], z0, z1)]:
            facet(f, (p[0], p[1], zlo), (q[0], q[1], zlo), (q[0], q[1], zhi))
            facet(f, (p[0], p[1], zlo), (q[0], q[1], zhi), (p[0], p[1], zhi))
        for z in (z0, z1):                                   # annular caps
            facet(f, (inn[i][0], inn[i][1], z), (out[i][0], out[i][1], z),
                  (out[j][0], out[j][1], z))
            facet(f, (inn[i][0], inn[i][1], z), (out[j][0], out[j][1], z),
                  (inn[j][0], inn[j][1], z))


def box(f, cx, cy, hx, hy, z0=0.0, z1=H):
    x0, x1, y0, y1 = cx - hx, cx + hx, cy - hy, cy + hy
    v = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
         (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    for a, b, c, d in [(0, 1, 2, 3), (7, 6, 5, 4), (0, 4, 5, 1),
                       (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]:
        facet(f, v[a], v[b], v[c])
        facet(f, v[a], v[c], v[d])


with open("tool.stl", "w") as f:
    f.write("solid winking_smiley\n")
    # solid back wall on top, so the tool is a real stamp (not a cookie cutter)
    cylinder(f, 0.0, 0.0, 28.0, H, H + 3.0)
    ring(f, 0.0, 0.0, 28.0, 23.0)               # outer ridge (the face border)
    cylinder(f, -10.5, 8.0, 4.5)                # left eye (open)
    box(f, 10.5, 8.5, 5.0, 1.9)                 # right eye (a wink)
    cx, cy, rad = 0.0, 3.0, 14.0                # smile: a downward arc of dots
    for k in range(13):
        ang = math.radians(216.0 + k * (108.0 / 12.0))
        cylinder(f, cx + rad * math.cos(ang), cy + rad * math.sin(ang), 2.8, n=24)
    f.write("endsolid winking_smiley\n")

print("wrote tool.stl: winking smiley stamp with outer ridge")
