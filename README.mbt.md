# MoonStrMech - Structural Mechanics Calculation Library

**MoonStrMech** is a comprehensive structural mechanics calculation library written in **Moonbit**, providing functions for structural analysis, beam theory, and civil engineering calculations.

## Features

### Support Reactions
- Simply supported beams (point loads, uniform loads)
- Cantilever beams (point loads, uniform loads)
- Overhanging beams

### Internal Forces
- Shear force calculations
- Bending moment calculations
- Beam deflection analysis

### Structural Analysis
- Truss analysis (triangle trusses, joint analysis)
- Stability analysis (Euler critical loads)
- Dynamic analysis (natural frequencies, damping)
- Continuous beam analysis

### Mathematical Utilities
- Linear algebra (matrix operations, determinants)
- Vector operations (2D/3D calculations)
- Geometric functions (polygon areas, distances)
- Numerical methods (approximations, safe divisions)

## Usage Examples

### Support Reactions

```moonbit
test "simply supported beam - point load" {
  // P=10kN at mid-span, L=4m
  let (ra, rb) = simply_supported_reactions_point(10.0, 2.0, 4.0)
  assert_eq(ra, 5.0)
  assert_eq(rb, 5.0)
}

test "simply supported beam - uniform load" {
  // q=5kN/m, L=6m
  let (ra, rb) = simply_supported_reactions_uniform(5.0, 6.0)
  assert_eq(ra, 15.0)
  assert_eq(rb, 15.0)
}

test "cantilever reactions" {
  // P=8kN at end, L=3m
  let (v, m) = cantilever_reactions_point(8.0, 3.0)
  assert_eq(v, 8.0)
  assert_eq(m, 24.0)

  // q=4kN/m, L=5m
  let (v2, m2) = cantilever_reactions_uniform(4.0, 5.0)
  assert_eq(v2, 20.0)
  assert_eq(m2, 50.0)
}
```

### Internal Forces and Deflection

```moonbit
test "internal forces - simply supported" {
  // P=10kN at 2m from left, L=5m
  let v_left = simply_supported_shear_point(10.0, 2.0, 5.0, 1.0)
  assert_eq(v_left, 6.0)
  let v_right = simply_supported_shear_point(10.0, 2.0, 5.0, 3.0)
  assert_eq(v_right, -4.0)

  let m_at_load = simply_supported_moment_point(10.0, 2.0, 5.0, 2.0)
  assert_eq(m_at_load, 12.0)
}

test "section properties" {
  // Rectangular section: b=0.3m, h=0.5m
  let i_rect = rect_inertia(0.3, 0.5)
  assert_eq(approx_equal(i_rect, 0.003125, 0.000001), true)

  let w_rect = rect_section_modulus(0.3, 0.5)
  assert_eq(approx_equal(w_rect, 0.0125, 0.000001), true)
}

test "deflection calculations" {
  // Simply supported beam, point load at mid-span
  // δ = PL³/(48EI)
  let d1 = simply_supported_deflection_point_mid(10.0, 4.0, 200.0, 0.001)
  assert_eq(approx_equal(d1, 66.667, 0.01), true)
}
```

### Stability and Dynamics

```moonbit
test "euler critical force" {
  // E=200GPa, I=0.0001m⁴, L=5m, μ=1 (both ends hinged)
  // Pcr = π²EI/(μL)² = 7.895 MN
  let pcr = euler_critical_force(200.0e9, 0.0001, 5.0, 1.0)
  assert_eq(approx_equal(pcr / 1.0e6, 7.895, 0.01), true)
}

test "stability factor - Q235 steel" {
  let phi = stability_factor_steel_q235(50.0)
  assert_eq(phi > 0.0, true)
}

test "dynamics" {
  // Natural frequency: k=1000N/m, m=10kg
  let omega = natural_frequency(1000.0, 10.0)
  assert_eq(omega, 10.0)

  // Damping ratio: c=40, k=1000, m=10
  let xi = damping_ratio(40.0, 1000.0, 10.0)
  assert_eq(approx_equal(xi, 0.2, 0.0001), true)
}
```

### Structural Analysis

```moonbit
test "truss analysis" {
  // Triangle truss: L=4m, h=3m, P=10kN
  let (n_left, n_right, n_bottom) = truss_triangle(4.0, 3.0, 10.0)
  assert_eq(approx_equal(n_left, -6.010, 0.01), true)
  assert_eq(approx_equal(n_right, -6.010, 0.01), true)
  assert_eq(approx_equal(n_bottom, 3.336, 0.01), true)
}

test "rotational stiffness" {
  // Rotational stiffness (far end fixed)
  let s_fixed = rotational_stiffness(200.0, 0.001, 4.0, true)
  assert_eq(s_fixed, 0.2)

  // Rotational stiffness (far end hinged)
  let s_hinged = rotational_stiffness(200.0, 0.001, 4.0, false)
  assert_eq(s_hinged, 0.15)
}

test "three moment equation" {
  let coeffs = three_moment_coefficients(4.0, 6.0)
  assert_eq(coeffs.0, 4.0)
}
```

### Mathematical Utilities

```moonbit
test "basic math tools" {
  assert_eq(square(3.0), 9.0)
  assert_eq(cube(2.0), 8.0)
  assert_eq(abs(-5.0), 5.0)
  assert_eq(max(3.0, 5.0), 5.0)
  assert_eq(min(3.0, 5.0), 3.0)
}

test "vector operations" {
  // Dot product
  assert_eq(dot2d(1.0, 0.0, 0.0, 1.0), 0.0)
  assert_eq(dot2d(1.0, 2.0, 3.0, 4.0), 11.0)

  // Cross product (2D)
  assert_eq(cross2d(1.0, 0.0, 0.0, 1.0), 1.0)

  // Distance
  assert_eq(distance2d(0.0, 0.0, 3.0, 4.0), 5.0)
}

test "determinants" {
  // 2x2 determinant
  assert_eq(det2x2(1.0, 2.0, 3.0, 4.0), -2.0)

  // 3x3 determinant
  assert_eq(det3x3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), 1.0)
}

test "solve linear systems" {
  // 2x2 system: x + 2y = 5, 3x + 4y = 11
  let (x, y) = solve2x2(1.0, 2.0, 3.0, 4.0, 5.0, 11.0)
  assert_eq(approx_equal(x, 1.0, 0.0001), true)
  assert_eq(approx_equal(y, 2.0, 0.0001), true)
}

test "polygon area" {
  // Square 2x2
  let sq = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
  assert_eq(polygon_area(sq), 4.0)

  // Triangle
  let tri = [(0.0, 0.0), (4.0, 0.0), (2.0, 3.0)]
  assert_eq(polygon_area(tri), 6.0)
}
```

### Force Method and Displacement Method

```moonbit
test "force method" {
  // Single redundant: δ₁p = 5, δ₁₁ = 2
  let x1 = force_method_1_redundant(5.0, 2.0)
  assert_eq(x1, -2.5)

  // Two redundants
  let (x1, x2) = force_method_2_redundant(1.0, 0.0, 1.0, 5.0, 3.0)
  assert_eq(approx_equal(x1, -5.0, 0.001), true)
  assert_eq(approx_equal(x2, -3.0, 0.001), true)
}

test "displacement method" {
  // Rotational stiffness (far end fixed)
  let s_fixed = rotational_stiffness(200.0, 0.001, 4.0, true)
  assert_eq(s_fixed, 0.2)

  // Linear stiffness
  let s_linear = linear_stiffness(200.0, 0.001, 4.0)
  assert_eq(s_linear, 0.05)
}
```

### Influence Lines

```moonbit
test "influence lines" {
  // Support reactions
  let eta_ra = influence_line_ra(2.0, 5.0)
  assert_eq(eta_ra, 0.6)

  let eta_rb = influence_line_rb(2.0, 5.0)
  assert_eq(eta_rb, 0.4)

  // Moment influence line (section at mid-span)
  let eta_m = influence_line_moment(4.0, 4.0, 10.0)
  assert_eq(eta_m, 2.4)
}
```

## Parameter Ranges

### Valid Input Ranges
- **Forces**: 0–10⁶ N (structural loads)
- **Lengths**: 0.1–100 m (beam spans and dimensions)
- **Elastic Modulus**: 10⁷–10¹¹ Pa (material stiffness)
- **Moments of Inertia**: 10⁻⁸–10⁻² m⁴ (section properties)

### Typical Engineering Values
- **Steel E**: 200 GPa (Young's modulus)
- **Concrete E**: 20–40 GPa (Young's modulus)
- **Wood E**: 8–15 GPa (Young's modulus)
- **Standard Beam Loads**: 1–100 kN/m (distributed loads)

## Testing

The project includes a comprehensive test suite covering all major functionalities:

```bash
moon test
```

### Test Coverage
- Support reactions for various beam types
- Internal force calculations (shear, moment)
- Beam deflection analysis
- Structural stability (Euler buckling)
- Dynamic analysis (natural frequency, damping)
- Mathematical utilities and geometric calculations
- Truss and arch structure analysis
- Force method and displacement method
- Influence line calculations

## Technical Details

### Engineering Standards
- **Structural Analysis**: Classical beam theory and structural mechanics
- **Euler Buckling**: Critical load analysis for slender members
- **Beam Deflection**: Integration methods for deflection calculations
- **SI Units**: Consistent use of International System of Units

### Applications
- **Civil Engineering**: Building structures, bridges, and foundations
- **Structural Design**: Beam analysis and design calculations
- **Educational Use**: Structural mechanics teaching and learning
- **Research Applications**: Structural analysis and simulation

## Notes

1. **Units**: All calculations use SI units (meters, kilograms, seconds, Newtons)
2. **Load Conventions**: Positive shear and moment follow standard conventions
3. **Boundary Conditions**: Standard support conditions assumed unless specified
4. **Material Properties**: Linear elastic behavior assumed
5. **Geometric Assumptions**: Euler-Bernoulli beam theory assumptions
6. **Validation**: Results should be verified against established structural analysis methods

## Version Information

The current version (0.1.0) implements **core structural mechanics calculation functions** including:

* Support reaction calculations for various beam types
* Internal force analysis (shear forces, bending moments)
* Beam deflection calculations using classical methods
* Structural stability analysis (Euler critical loads)
* Dynamic analysis (natural frequencies, damping)
* Mathematical utilities for structural calculations
* Truss and arch structure analysis
* Force method and displacement method implementations
* Influence line calculations for moving loads

The library is actively developed and aims to provide comprehensive coverage of structural mechanics while being optimized for MoonBit.
