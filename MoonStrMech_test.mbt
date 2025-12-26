// MoonStrMech - 测试代码

///|
test "basic math tools" {
  // square, cube, pow4
  assert_eq(square(3.0), 9.0)
  assert_eq(cube(2.0), 8.0)
  assert_eq(pow4(2.0), 16.0)
  // abs
  assert_eq(abs(-5.0), 5.0)
  assert_eq(abs(5.0), 5.0)
  // max, min
  assert_eq(max(3.0, 5.0), 5.0)
  assert_eq(min(3.0, 5.0), 3.0)
}

///|
test "safe_div and approx_equal" {
  assert_eq(safe_div(10.0, 2.0, 0.0), 5.0)
  assert_eq(safe_div(10.0, 0.0, -1.0), -1.0)
  assert_eq(approx_equal(1.0, 1.0000000001, 0.001), true)
  assert_eq(approx_equal(1.0, 2.0, 0.1), false)
}

///|
test "vector operations" {
  // dot product
  assert_eq(dot2d(1.0, 0.0, 0.0, 1.0), 0.0)
  assert_eq(dot2d(1.0, 2.0, 3.0, 4.0), 11.0)
  // cross product (2D returns scalar)
  assert_eq(cross2d(1.0, 0.0, 0.0, 1.0), 1.0)
  // distance
  assert_eq(distance2d(0.0, 0.0, 3.0, 4.0), 5.0)
}

///|
test "determinants" {
  // 2x2 det
  assert_eq(det2x2(1.0, 2.0, 3.0, 4.0), -2.0)
  // 3x3 det (identity has det = 1)
  assert_eq(det3x3(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), 1.0)
}

///|
test "solve linear systems" {
  // 2x2: x + 2y = 5, 3x + 4y = 11 => x = 1, y = 2
  let (x, y) = solve2x2(1.0, 2.0, 3.0, 4.0, 5.0, 11.0)
  assert_eq(approx_equal(x, 1.0, 0.0001), true)
  assert_eq(approx_equal(y, 2.0, 0.0001), true)
}

///|
test "polygon area" {
  // square 2x2
  let sq = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
  assert_eq(polygon_area(sq), 4.0)
  // triangle
  let tri = [(0.0, 0.0), (4.0, 0.0), (2.0, 3.0)]
  assert_eq(polygon_area(tri), 6.0)
}

///|
test "simply supported beam reactions - point load" {
  // P=10kN at mid-span, L=4m
  let (ra, rb) = simply_supported_reactions_point(10.0, 2.0, 4.0)
  assert_eq(ra, 5.0)
  assert_eq(rb, 5.0)
  // P=12kN at 1m from left, L=4m
  let (ra2, rb2) = simply_supported_reactions_point(12.0, 1.0, 4.0)
  assert_eq(ra2, 9.0)
  assert_eq(rb2, 3.0)
}

///|
test "simply supported beam reactions - uniform load" {
  // q=5kN/m, L=6m
  let (ra, rb) = simply_supported_reactions_uniform(5.0, 6.0)
  assert_eq(ra, 15.0)
  assert_eq(rb, 15.0)
}

///|
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

///|
test "internal forces - simply supported point load" {
  // P=10kN at 2m from left, L=5m
  // RA = 10*3/5 = 6kN
  let v_left = simply_supported_shear_point(10.0, 2.0, 5.0, 1.0)
  assert_eq(v_left, 6.0)
  let v_right = simply_supported_shear_point(10.0, 2.0, 5.0, 3.0)
  assert_eq(v_right, -4.0)
  // moment at load point
  let m_at_load = simply_supported_moment_point(10.0, 2.0, 5.0, 2.0)
  assert_eq(m_at_load, 12.0)
}

///|
test "internal forces - uniform load" {
  // q=8kN/m, L=4m, max moment at mid-span
  let m_max = simply_supported_max_moment_uniform(8.0, 4.0)
  assert_eq(m_max, 16.0)
  // shear at left support
  let v_0 = simply_supported_shear_uniform(8.0, 4.0, 0.0)
  assert_eq(v_0, 16.0)
  // shear at mid-span
  let v_mid = simply_supported_shear_uniform(8.0, 4.0, 2.0)
  assert_eq(v_mid, 0.0)
}

///|
test "section properties" {
  // rectangular section 0.3m x 0.5m
  let i_rect = rect_inertia(0.3, 0.5)
  assert_eq(approx_equal(i_rect, 0.003125, 0.000001), true)
  // section modulus
  let w_rect = rect_section_modulus(0.3, 0.5)
  assert_eq(approx_equal(w_rect, 0.0125, 0.000001), true)
}

///|
test "deflection calculations" {
  // Simply supported, point load at mid
  // δ = PL³/(48EI)
  // P=10, L=4, E=200, I=0.001 => 10*64/(48*200*0.001) = 640/9.6 = 66.67
  let d1 = simply_supported_deflection_point_mid(10.0, 4.0, 200.0, 0.001)
  assert_eq(approx_equal(d1, 66.667, 0.01), true)
  // Cantilever, point load at end
  // δ = PL³/(3EI)
  // P=5, L=2, E=100, I=0.002 => 5*8/(3*100*0.002) = 40/0.6 = 66.67
  let d2 = cantilever_deflection_point_end(5.0, 2.0, 100.0, 0.002)
  assert_eq(approx_equal(d2, 66.667, 0.01), true)
}

///|
test "truss analysis" {
  // Simple triangle truss: L=4m, h=3m, P=10kN
  // slant length = sqrt((L/2)^2 + h^2) = sqrt(4 + 9) = sqrt(13) ≈ 3.606m
  // sin = 3/3.606 ≈ 0.832, cos = 2/3.606 ≈ 0.555
  // N_slant = -P/(2*sin) = -10/(2*0.832) ≈ -6.010
  // N_bottom = -N_slant * cos = 6.010 * 0.555 ≈ 3.336
  let (n_left, n_right, n_bottom) = truss_triangle(4.0, 3.0, 10.0)
  assert_eq(approx_equal(n_left, -6.010, 0.01), true)
  assert_eq(approx_equal(n_right, -6.010, 0.01), true)
  assert_eq(approx_equal(n_bottom, 3.336, 0.01), true)
}

///|
test "rotational stiffness and distribution" {
  // Linear stiffness
  let ic = linear_stiffness(200.0, 0.001, 4.0)
  assert_eq(ic, 0.05)
  // Rotational stiffness (far end fixed)
  let s_fixed = rotational_stiffness(200.0, 0.001, 4.0, true)
  assert_eq(s_fixed, 0.2)
  // Rotational stiffness (far end hinged)
  let s_hinged = rotational_stiffness(200.0, 0.001, 4.0, false)
  assert_eq(s_hinged, 0.15)
  // Distribution factors
  let factors = distribution_factors([0.2, 0.3, 0.5])
  assert_eq(approx_equal(factors[0], 0.2, 0.0001), true)
  assert_eq(approx_equal(factors[1], 0.3, 0.0001), true)
  assert_eq(approx_equal(factors[2], 0.5, 0.0001), true)
}

///|
test "fixed end moments" {
  // Uniform load, both ends fixed
  // M = -qL²/12
  let (mab, mba) = fixed_end_moment_uniform(12.0, 6.0)
  assert_eq(mab, -36.0)
  assert_eq(mba, -36.0)
  // Point load at mid-span
  let (mab2, mba2) = fixed_end_moment_point_mid(24.0, 8.0)
  assert_eq(mab2, -24.0)
  assert_eq(mba2, 24.0)
}

///|
test "influence lines" {
  // RA influence at x=0 should be 1, at x=L should be 0
  assert_eq(influence_line_ra(0.0, 10.0), 1.0)
  assert_eq(influence_line_ra(10.0, 10.0), 0.0)
  assert_eq(influence_line_ra(5.0, 10.0), 0.5)
  // Moment influence at section a=4, L=10
  // at x=4: η = 4*6/10 = 2.4
  let eta_m = influence_line_moment(4.0, 4.0, 10.0)
  assert_eq(eta_m, 2.4)
}

///|
test "euler critical force" {
  // E=200GPa, I=0.0001m⁴, L=5m, μ=1 (both ends hinged)
  // Pcr = π²EI/(μL)² = π²*200*0.0001/25 = 0.0079 * 10^9 = 7.895 * 10^6
  let pcr = euler_critical_force(200.0e9, 0.0001, 5.0, 1.0)
  assert_eq(approx_equal(pcr / 1.0e6, 7.895, 0.01), true)
}

///|
test "stability factor" {
  // λ = 50
  let phi1 = stability_factor_steel_q235(50.0)
  assert_eq(approx_equal(phi1, 0.666, 0.001), true)
  // λ = 100
  let phi2 = stability_factor_steel_q235(100.0)
  assert_eq(approx_equal(phi2, 0.332, 0.001), true)
}

///|
test "arch structures" {
  // Three-hinged arch: q=10kN/m, L=20m, f=5m
  // H = qL²/(8f) = 10*400/40 = 100kN
  let h = three_hinged_arch_thrust_uniform(10.0, 20.0, 5.0)
  assert_eq(h, 100.0)
  // Parabolic arch y at mid-span
  // y(L/2) = 4f*(L/2)*(L/2)/L² = f
  let y_mid = parabolic_arch_y(10.0, 20.0, 5.0)
  assert_eq(y_mid, 5.0)
}

///|
test "strain energy" {
  // Bending: M=100kNm, L=4m, E=200GPa, I=0.001m⁴
  // U = M²L/(2EI) = 10000*4/(2*200e9*0.001) = 40000/(4e8) = 0.0001 kNm = 0.1 J (if consistent units)
  let u_bend = bending_strain_energy(100.0, 4.0, 200.0e9, 0.001)
  assert_eq(approx_equal(u_bend, 0.0001, 0.00001), true)
}

///|
test "dynamics" {
  // Natural frequency: k=1000N/m, m=10kg
  // ω = sqrt(k/m) = sqrt(100) = 10 rad/s
  let omega = natural_frequency(1000.0, 10.0)
  assert_eq(omega, 10.0)
  // Natural period T = 2π/ω
  let t = natural_period(1000.0, 10.0)
  assert_eq(approx_equal(t, 0.6283, 0.0001), true)
  // Damping ratio
  let xi = damping_ratio(40.0, 1000.0, 10.0)
  assert_eq(approx_equal(xi, 0.2, 0.0001), true)
}

///|
test "continuous beam" {
  // Two-span continuous beam, equal spans, uniform load
  // MB = -qL²/8
  let mb = continuous_beam_2span_uniform_mid_moment(10.0, 6.0)
  assert_eq(mb, -45.0)
}

///|
test "graph multiplication" {
  // Two triangles with same base
  let result = graph_multiply_triangles(6.0, 4.0, 3.0)
  assert_eq(result, 24.0)
  // Triangle and rectangle
  let result2 = graph_multiply_tri_rect(6.0, 2.0, 4.0)
  assert_eq(result2, 24.0)
}

///|
test "force method" {
  // One redundant: delta_1p = -10, delta_11 = 2
  // X1 = -(-10)/2 = 5
  let x1 = force_method_1_redundant(-10.0, 2.0)
  assert_eq(x1, 5.0)
}

///|
test "slenderness and radius of gyration" {
  // I=0.0001, A=0.01, L=5, μ=1
  // i = sqrt(I/A) = sqrt(0.01) = 0.1
  // λ = μL/i = 1*5/0.1 = 50
  let i = radius_of_gyration(0.0001, 0.01)
  assert_eq(i, 0.1)
  let lambda = slenderness_ratio(5.0, 1.0, 0.0001, 0.01)
  assert_eq(lambda, 50.0)
}

///|
test "virtual work" {
  // Bending contribution: M_actual=10, M_virtual=1, L=4, E=200, I=0.001
  // Δ = M*M̄*L/(EI) = 10*1*4/(200*0.001) = 40/0.2 = 200
  let delta = virtual_work_bending(10.0, 1.0, 4.0, 200.0, 0.001)
  assert_eq(delta, 200.0)
}
