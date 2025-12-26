# MoonStrMech - Structural Mechanics Library

A comprehensive structural mechanics calculation library in Moonbit for civil engineering applications.

## Features

- **Beam Analysis**: Support reactions, shear forces, bending moments, deflections
- **Structural Analysis**: Truss analysis, stability, dynamics
- **Mathematical Tools**: Linear algebra, vectors, geometry functions

## Quick Start

```moonbit
// Support reactions
let (ra, rb) = simply_supported_reactions_point(10.0, 2.0, 4.0)

// Internal forces
let shear = simply_supported_shear_point(10.0, 2.0, 5.0, 1.0)
let moment = simply_supported_moment_point(10.0, 2.0, 5.0, 2.0)

// Stability analysis
let critical_force = euler_critical_force(200.0e9, 0.0001, 5.0, 1.0)

// Truss analysis
let (n1, n2, n3) = truss_triangle(4.0, 3.0, 10.0)
```

## Testing

```bash
moon test
```

## Technical Details

- **Units**: SI (N, m, Pa)
- **Theory**: Classical beam theory, Euler buckling
- **Applications**: Structural design, civil engineering, education

## Version

v0.1.0 - Core structural mechanics functions for Moonbit
