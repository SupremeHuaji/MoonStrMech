# MoonStrMech - Structural Mechanics Calculation Library

**MoonStrMech** is a comprehensive structural mechanics calculation library written in **Moonbit**, providing functions for structural analysis, beam theory, and civil engineering calculations.

---

## Features

### Support Reactions

* **Simply Supported Beams**: Reactions for point loads, uniform loads, and triangular loads
* **Cantilever Beams**: Fixed-end reactions for point loads and uniform loads
* **Overhanging Beams**: Support reactions for complex beam configurations

---

### Internal Forces

* **Shear Forces**: Calculation of shear force diagrams for various load conditions
* **Bending Moments**: Bending moment calculations and diagrams
* **Deflections**: Beam deflection analysis and calculations

---

### Structural Analysis

* **Truss Analysis**: Force calculations in truss structures
* **Arch Structures**: Three-hinged arch analysis
* **Continuous Beams**: Multi-span beam analysis
* **Stability Analysis**: Euler critical loads and stability factors

---

### Mathematical Utilities

* **Linear Algebra**: Matrix operations, determinants, and system solving
* **Vector Operations**: 2D/3D vector calculations, dot products, cross products
* **Geometric Functions**: Polygon area calculations, distance computations
* **Numerical Methods**: Approximation functions and safe divisions

---

## Usage

### Support Reactions

```moonbit
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

test "simply supported beam reactions - uniform load" {
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

---

### Internal Forces and Deflection

```moonbit
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

test "section properties" {
  // rectangular section 0.3m x 0.5m
  let i_rect = rect_inertia(0.3, 0.5)
  assert_eq(approx_equal(i_rect, 0.003125, 0.000001), true)
  // section modulus
  let w_rect = rect_section_modulus(0.3, 0.5)
  assert_eq(approx_equal(w_rect, 0.0125, 0.000001), true)
}

test "deflection calculations" {
  // Simply supported, point load at mid
  // δ = PL³/(48EI)
  let d1 = simply_supported_deflection_point_mid(10.0, 4.0, 200.0, 0.001)
  assert_eq(approx_equal(d1, 66.667, 0.01), true)
}
```

---

### Stability and Dynamics

```moonbit
test "euler critical force" {
  // E=200GPa, I=0.0001m⁴, L=5m, μ=1 (both ends hinged)
  // Pcr = π²EI/(μL)² = 7.895 MN
  let pcr = euler_critical_force(200.0e9, 0.0001, 5.0, 1.0)
  assert_eq(approx_equal(pcr / 1.0e6, 7.895, 0.01), true)
}

test "dynamics" {
  // Natural frequency: k=1000N/m, m=10kg
  // ω = sqrt(k/m) = 10 rad/s
  let omega = natural_frequency(1000.0, 10.0)
  assert_eq(omega, 10.0)
  // Damping ratio
  let xi = damping_ratio(40.0, 1000.0, 10.0)
  assert_eq(approx_equal(xi, 0.2, 0.0001), true)
}
```

---

### Structural Analysis

```moonbit
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
```

---

### Mathematical Utilities

```moonbit
test "basic math tools" {
  assert_eq(square(3.0), 9.0)
  assert_eq(cube(2.0), 8.0)
  assert_eq(abs(-5.0), 5.0)
  assert_eq(max(3.0, 5.0), 5.0)
}

test "vector operations" {
  assert_eq(dot2d(1.0, 0.0, 0.0, 1.0), 0.0)
  assert_eq(cross2d(1.0, 0.0, 0.0, 1.0), 1.0)
  assert_eq(distance2d(0.0, 0.0, 3.0, 4.0), 5.0)
}

test "solve linear systems" {
  // 2x2: x + 2y = 5, 3x + 4y = 11 => x = 1, y = 2
  let (x, y) = solve2x2(1.0, 2.0, 3.0, 4.0, 5.0, 11.0)
  assert_eq(approx_equal(x, 1.0, 0.0001), true)
  assert_eq(approx_equal(y, 2.0, 0.0001), true)
}

test "polygon area" {
  let sq = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
  assert_eq(polygon_area(sq), 4.0)
}
```

---

## Parameter Ranges

### Valid Input Ranges

* **Forces**: 0–10⁶ N (structural loads)
* **Lengths**: 0.1–100 m (beam spans and dimensions)
* **Elastic Modulus**: 10⁷–10¹¹ Pa (material stiffness)
* **Moments of Inertia**: 10⁻⁸–10⁻² m⁴ (section properties)

---

### Typical Engineering Values

* **Steel E**: 200 GPa (Young's modulus)
* **Concrete E**: 20–40 GPa (Young's modulus)
* **Wood E**: 8–15 GPa (Young's modulus)
* **Standard Beam Loads**: 1–100 kN/m (distributed loads)

---

## Testing

The project includes a comprehensive test suite covering all major functionalities:

```bash
moon test
```

### Test Coverage

* Support reactions for various beam types
* Internal force calculations (shear, moment)
* Beam deflection analysis
* Structural stability (Euler buckling)
* Dynamic analysis (natural frequency, damping)
* Mathematical utilities and geometric calculations
* Truss and arch structure analysis

---

## Technical Details

### Engineering Standards

* **Structural Analysis**: Classical beam theory and structural mechanics
* **Euler Buckling**: Critical load analysis for slender members
* **Beam Deflection**: Integration methods for deflection calculations
* **SI Units**: Consistent use of International System of Units

---

### Applications

* **Civil Engineering**: Building structures, bridges, and foundations
* **Structural Design**: Beam analysis and design calculations
* **Educational Use**: Structural mechanics teaching and learning
* **Research Applications**: Structural analysis and simulation

---

## Notes

1. **Units**: All calculations use SI units (meters, kilograms, seconds, Newtons)
2. **Load Conventions**: Positive shear and moment follow standard conventions
3. **Boundary Conditions**: Standard support conditions assumed unless specified
4. **Material Properties**: Linear elastic behavior assumed
5. **Geometric Assumptions**: Euler-Bernoulli beam theory assumptions
6. **Validation**: Results should be verified against established structural analysis methods

---

## Version Information

The current version (0.1.0) implements **core structural mechanics calculation functions** including:

* Support reaction calculations for various beam types
* Internal force analysis (shear forces, bending moments)
* Beam deflection calculations using classical methods
* Structural stability analysis (Euler critical loads)
* Dynamic analysis (natural frequencies, damping)
* Mathematical utilities for structural calculations
* Truss and arch structure analysis

The library is actively developed and aims to provide comprehensive coverage of structural mechanics while being optimized for MoonBit.
  l : Double,
  x : Double,
) -> Double {
  if l <= 0.0 {
    return 0.0
  }
  let (ra, _) = simply_supported_reactions_point(p, a, l)
  if x < a {
    ra
  } else {
    ra - p
  }
}


///|
/// 简支梁集中力作用下任意截面的弯矩
/// P 作用于距左端 a 处，求距左端 x 处的弯矩
pub fn simply_supported_moment_point(
  p : Double,
  a : Double,
  l : Double,
  x : Double,
) -> Double {
  if l <= 0.0 {
    return 0.0
  }
  let (ra, _) = simply_supported_reactions_point(p, a, l)
  if x <= a {
    ra * x
  } else {
    ra * x - p * (x - a)
  }
}

///|
/// 简支梁均布荷载作用下任意截面的剪力
pub fn simply_supported_shear_uniform(
  q : Double,
  l : Double,
  x : Double,
) -> Double {
  if l <= 0.0 {
    return 0.0
  }
  let (ra, _) = simply_supported_reactions_uniform(q, l)
  ra - q * x
}

///|
/// 简支梁均布荷载作用下任意截面的弯矩
pub fn simply_supported_moment_uniform(
  q : Double,
  l : Double,
  x : Double,
) -> Double {
  if l <= 0.0 {
    return 0.0
  }
  let (ra, _) = simply_supported_reactions_uniform(q, l)
  ra * x - q * x * x / 2.0
}

///|
/// 简支梁均布荷载下的最大弯矩 (跨中)
pub fn simply_supported_max_moment_uniform(q : Double, l : Double) -> Double {
  q * l * l / 8.0
}

///|
/// 悬臂梁端部集中力作用下任意截面的剪力
/// x 从固定端量起
pub fn cantilever_shear_point(p : Double, l : Double, x : Double) -> Double {
  if x < l {
    p
  } else {
    0.0
  }
}

///|
/// 悬臂梁端部集中力作用下任意截面的弯矩
/// x 从固定端量起，弯矩为负(上部受拉)
pub fn cantilever_moment_point(p : Double, l : Double, x : Double) -> Double {
  if x <= l {
    -p * (l - x)
  } else {
    0.0
  }
}

///|
/// 悬臂梁均布荷载作用下任意截面的剪力
pub fn cantilever_shear_uniform(q : Double, l : Double, x : Double) -> Double {
  if x <= l {
    q * (l - x)
  } else {
    0.0
  }
}

///|
/// 悬臂梁均布荷载作用下任意截面的弯矩
pub fn cantilever_moment_uniform(q : Double, l : Double, x : Double) -> Double {
  if x <= l {
    let remaining = l - x
    -q * remaining * remaining / 2.0
  } else {
    0.0
  }
}

// ==================== 截面特性 ====================

///|
/// 矩形截面惯性矩 Ix = bh³/12
pub fn rect_inertia(b : Double, h : Double) -> Double {
  b * cube(h) / 12.0
}

///|
/// 矩形截面抗弯模量 Wx = bh²/6
pub fn rect_section_modulus(b : Double, h : Double) -> Double {
  b * square(h) / 6.0
}

///|
/// 圆形截面惯性矩 I = πr⁴/4
pub fn circle_inertia(r : Double) -> Double {
  PI * pow4(r) / 4.0
}

///|
/// 圆形截面抗弯模量 W = πr³/4
pub fn circle_section_modulus(r : Double) -> Double {
  PI * cube(r) / 4.0
}

///|
/// 工字形截面惯性矩 (简化)
/// b - 翼缘宽度, h - 总高度, tf - 翼缘厚度, tw - 腹板厚度
pub fn i_section_inertia(
  b : Double,
  h : Double,
  tf : Double,
  tw : Double,
) -> Double {
  // Ix = (b*h³ - (b-tw)*(h-2tf)³) / 12
  let outer = b * cube(h) / 12.0
  let inner_h = h - 2.0 * tf
  let inner_b = b - tw
  let inner = if inner_h > 0.0 && inner_b > 0.0 {
    inner_b * cube(inner_h) / 12.0
  } else {
    0.0
  }
  outer - inner
}

///|
/// 空心矩形截面惯性矩
/// B, H - 外尺寸，b, h - 内空尺寸
pub fn hollow_rect_inertia(
  b_outer : Double,
  h_outer : Double,
  b_inner : Double,
  h_inner : Double,
) -> Double {
  rect_inertia(b_outer, h_outer) - rect_inertia(b_inner, h_inner)
}

///|
/// 空心圆形截面惯性矩
pub fn hollow_circle_inertia(r_outer : Double, r_inner : Double) -> Double {
  PI * (pow4(r_outer) - pow4(r_inner)) / 4.0
}

///|
/// 回转半径 i = sqrt(I/A)
pub fn radius_of_gyration(inertia : Double, area : Double) -> Double {
  if area <= 0.0 {
    0.0
  } else {
    (inertia / area).sqrt()
  }
}

///|
/// 长细比 λ = μL/i
pub fn slenderness_ratio(
  length : Double,
  mu : Double,
  inertia : Double,
  area : Double,
) -> Double {
  let i = radius_of_gyration(inertia, area)
  if i <= 0.0 {
    0.0
  } else {
    mu * length / i
  }
}

// ==================== 位移计算 ====================

///|
/// 简支梁中点集中力作用下的跨中挠度
/// δ = PL³/(48EI)
pub fn simply_supported_deflection_point_mid(
  p : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  p * cube(l) / (48.0 * e * i)
}

///|
/// 简支梁均布荷载作用下的跨中挠度
/// δ = 5qL⁴/(384EI)
pub fn simply_supported_deflection_uniform_mid(
  q : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  5.0 * q * pow4(l) / (384.0 * e * i)
}

///|
/// 悬臂梁端部集中力作用下的自由端挠度
/// δ = PL³/(3EI)
pub fn cantilever_deflection_point_end(
  p : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  p * cube(l) / (3.0 * e * i)
}

///|
/// 悬臂梁均布荷载作用下的自由端挠度
/// δ = qL⁴/(8EI)
pub fn cantilever_deflection_uniform_end(
  q : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  q * pow4(l) / (8.0 * e * i)
}

///|
/// 悬臂梁端部集中力作用下的自由端转角
/// θ = PL²/(2EI)
pub fn cantilever_slope_point_end(
  p : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  p * square(l) / (2.0 * e * i)
}

///|
/// 简支梁端部转角 (均布荷载)
/// θ = qL³/(24EI)
pub fn simply_supported_slope_uniform_end(
  q : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  q * cube(l) / (24.0 * e * i)
}

// ==================== 桁架分析 ====================

///|
/// 计算桁架节点法 - 两杆交汇节点
/// 给定两杆角度和外力，求两杆内力
/// 角度相对于水平方向，逆时针为正
/// 返回 (杆1内力, 杆2内力)，正为拉力
pub fn truss_joint_2bars(
  fx : Double,
  fy : Double,
  angle1_deg : Double,
  angle2_deg : Double,
) -> (Double, Double) {
  let a1 = deg_to_rad(angle1_deg)
  let a2 = deg_to_rad(angle2_deg)
  let cos1 = cos_approx(a1)
  let sin1 = sin_approx(a1)
  let cos2 = cos_approx(a2)
  let sin2 = sin_approx(a2)
  // 平衡方程:
  // N1*cos1 + N2*cos2 + Fx = 0
  // N1*sin1 + N2*sin2 + Fy = 0
  let (n1, n2) = solve2x2(cos1, cos2, sin1, sin2, -fx, -fy)
  (n1, n2)
}

///|
/// 简单三角形桁架分析
/// 三节点桁架，底边水平，顶点在上方
/// 给定底边长度 L，高度 h，顶点荷载 P (向下为正)
/// 返回 (左斜杆内力, 右斜杆内力, 底杆内力)
pub fn truss_triangle(
  l : Double,
  h : Double,
  p : Double,
) -> (Double, Double, Double) {
  if l <= 0.0 || h <= 0.0 {
    return (0.0, 0.0, 0.0)
  }
  // 左斜杆和右斜杆的长度和角度
  let half_l = l / 2.0
  let slant_len = (half_l * half_l + h * h).sqrt()
  let sin_theta = h / slant_len
  let cos_theta = half_l / slant_len
  // 顶点平衡:
  // N_left * sin + N_right * sin = P
  // N_left * cos = N_right * cos (对称)
  // 所以 N_left = N_right = P / (2*sin)
  let n_slant = -p / (2.0 * sin_theta) // 压力为负
  // 底杆：取左节点平衡
  // 左支座反力 = P/2 (向上)
  // N_left * cos + N_bottom = 0
  let n_bottom = -n_slant * cos_theta // 拉力
  (n_slant, n_slant, n_bottom)
}

// ==================== 力法基础 ====================

///|
/// 计算单位荷载图乘法积分
/// 两个三角形分布的图乘 (共基边)
/// M1: 最大值 m1，M2: 最大值 m2，跨度 L
/// 结果 = (1/3) * L * m1 * m2
pub fn graph_multiply_triangles(m1 : Double, m2 : Double, l : Double) -> Double {
  l * m1 * m2 / 3.0
}

///|
/// 计算图乘 - 三角形与矩形
/// 三角形最大值 m_tri，矩形高度 m_rect，跨度 L
/// 结果 = (1/2) * L * m_tri * m_rect
pub fn graph_multiply_tri_rect(
  m_tri : Double,
  m_rect : Double,
  l : Double,
) -> Double {
  l * m_tri * m_rect / 2.0
}

///|
/// 计算图乘 - 两个矩形
/// 结果 = L * m1 * m2
pub fn graph_multiply_rects(m1 : Double, m2 : Double, l : Double) -> Double {
  l * m1 * m2
}

///|
/// 计算图乘 - 抛物线与三角形 (简支梁均布弯矩与单位荷载三角形)
/// 抛物线最大值 m_para (跨中)，三角形最大值 m_tri，跨度 L
/// 积分 = (1/3) * L * m_para * m_tri (当三角形顶点在跨中时)
/// 一般情况: (2/3) * L * m_para * m_tri_mid (m_tri_mid是三角形在跨中的值)
pub fn graph_multiply_para_tri(
  m_para : Double,
  m_tri : Double,
  l : Double,
) -> Double {
  2.0 * l * m_para * m_tri / 3.0
}

///|
/// 力法：一次超静定梁的多余力计算
/// delta_1p: 基本结构在荷载作用下沿多余约束方向的位移
/// delta_11: 基本结构在单位多余力作用下沿多余约束方向的位移
/// 返回多余力 X1
pub fn force_method_1_redundant(delta_1p : Double, delta_11 : Double) -> Double {
  if abs(delta_11) < EPSILON {
    0.0
  } else {
    -delta_1p / delta_11
  }
}

///|
/// 力法：二次超静定结构的多余力计算
/// 输入柔度矩阵元素和荷载位移
/// 返回 (X1, X2)
pub fn force_method_2_redundant(
  d11 : Double,
  d12 : Double,
  d22 : Double,
  d1p : Double,
  d2p : Double,
) -> (Double, Double) {
  // [d11 d12] [X1]   [-d1p]
  // [d12 d22] [X2] = [-d2p]
  let (x1, x2) = solve2x2(d11, d12, d12, d22, -d1p, -d2p)
  (x1, x2)
}

// ==================== 位移法基础 ====================

///|
/// 等截面直杆的转动刚度
/// 远端固定: S = 4EI/L
/// 远端铰接: S = 3EI/L
pub fn rotational_stiffness(
  e : Double,
  i : Double,
  l : Double,
  far_end_fixed : Bool,
) -> Double {
  if l <= 0.0 {
    return 0.0
  }
  if far_end_fixed {
    4.0 * e * i / l
  } else {
    3.0 * e * i / l
  }
}

///|
/// 等截面直杆的线刚度
/// i_c = EI/L
pub fn linear_stiffness(e : Double, i : Double, l : Double) -> Double {
  if l <= 0.0 {
    0.0
  } else {
    e * i / l
  }
}

///|
/// 传递系数
/// 远端固定: C = 0.5
/// 远端铰接: C = 0
pub fn carryover_factor(far_end_fixed : Bool) -> Double {
  if far_end_fixed {
    0.5
  } else {
    0.0
  }
}

///|
/// 固端弯矩 - 简支梁两端固定，中点集中力
/// M_AB = -PL/8, M_BA = PL/8
pub fn fixed_end_moment_point_mid(p : Double, l : Double) -> (Double, Double) {
  let m = p * l / 8.0
  (-m, m)
}

///|
/// 固端弯矩 - 两端固定，集中力作用于距A端 a 处
/// M_AB = -Pab²/L², M_BA = Pa²b/L²
pub fn fixed_end_moment_point(
  p : Double,
  a : Double,
  l : Double,
) -> (Double, Double) {
  if l <= 0.0 {
    return (0.0, 0.0)
  }
  let b = l - a
  let m_ab = -p * a * b * b / (l * l)
  let m_ba = p * a * a * b / (l * l)
  (m_ab, m_ba)
}

///|
/// 固端弯矩 - 两端固定，均布荷载
/// M_AB = M_BA = -qL²/12
pub fn fixed_end_moment_uniform(q : Double, l : Double) -> (Double, Double) {
  let m = -q * l * l / 12.0
  (m, m)
}

///|
/// 固端弯矩 - 两端固定，三角形荷载 (从A端0增到B端q_max)
/// M_AB = -qL²/30, M_BA = qL²/20
pub fn fixed_end_moment_triangular(
  q_max : Double,
  l : Double,
) -> (Double, Double) {
  let m_ab = -q_max * l * l / 30.0
  let m_ba = q_max * l * l / 20.0
  (m_ab, m_ba)
}

// ==================== 力矩分配法 ====================

///|
/// 计算分配系数
/// 给定节点各杆件的转动刚度，返回各杆件的分配系数
pub fn distribution_factors(stiffnesses : Array[Double]) -> Array[Double] {
  let total = array_sum(stiffnesses)
  if abs(total) < EPSILON {
    // 返回零数组
    let result : Array[Double] = []
    for i = 0; i < stiffnesses.length(); i = i + 1 {
      result.push(0.0)
    }
    return result
  }
  let result : Array[Double] = []
  for i = 0; i < stiffnesses.length(); i = i + 1 {
    result.push(stiffnesses[i] / total)
  }
  result
}

///|
/// 单节点力矩分配 - 一次分配
/// unbalanced_moment: 不平衡力矩
/// factors: 分配系数数组
/// 返回各杆件分配的力矩
pub fn distribute_moment(
  unbalanced_moment : Double,
  factors : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = []
  for i = 0; i < factors.length(); i = i + 1 {
    result.push(-unbalanced_moment * factors[i])
  }
  result
}

///|
/// 传递力矩计算
/// distributed: 分配的力矩
/// carryover: 传递系数
/// 返回传递到远端的力矩
pub fn carryover_moment(distributed : Double, carryover : Double) -> Double {
  distributed * carryover
}

// ==================== 影响线 ====================

///|
/// 简支梁支座反力影响线纵标
/// 求 RA 在荷载位于 x 处时的影响线值
/// η_RA(x) = (L - x) / L
pub fn influence_line_ra(x : Double, l : Double) -> Double {
  if l <= 0.0 || x < 0.0 || x > l {
    0.0
  } else {
    (l - x) / l
  }
}

///|
/// 简支梁支座反力影响线纵标
/// 求 RB 在荷载位于 x 处时的影响线值
/// η_RB(x) = x / L
pub fn influence_line_rb(x : Double, l : Double) -> Double {
  if l <= 0.0 || x < 0.0 || x > l {
    0.0
  } else {
    x / l
  }
}

///|
/// 简支梁某截面弯矩影响线纵标
/// 截面位于距左端 a 处，荷载位于 x 处
pub fn influence_line_moment(x : Double, a : Double, l : Double) -> Double {
  if l <= 0.0 || x < 0.0 || x > l || a < 0.0 || a > l {
    return 0.0
  }
  let b = l - a
  if x <= a {
    x * b / l
  } else {
    a * (l - x) / l
  }
}

///|
/// 简支梁某截面剪力影响线纵标
/// 截面位于距左端 a 处，荷载位于 x 处
pub fn influence_line_shear(x : Double, a : Double, l : Double) -> Double {
  if l <= 0.0 || x < 0.0 || x > l || a < 0.0 || a > l {
    return 0.0
  }
  if x < a {
    -x / l // 左侧
  } else {
    (l - x) / l // 右侧
  }
}

///|
/// 利用影响线计算均布荷载作用效应
/// 影响线面积法: S = q * Ω
/// omega: 影响线包围面积
pub fn effect_from_influence_area(q : Double, omega : Double) -> Double {
  q * omega
}

///|
/// 三角形影响线面积
pub fn influence_area_triangle(base : Double, height : Double) -> Double {
  base * height / 2.0
}

// ==================== 稳定性计算 ====================

///|
/// 欧拉临界力
/// P_cr = π²EI / (μL)²
pub fn euler_critical_force(
  e : Double,
  i : Double,
  l : Double,
  mu : Double,
) -> Double {
  if l <= 0.0 || mu <= 0.0 {
    return 0.0
  }
  PI * PI * e * i / square(mu * l)
}

///|
/// 欧拉临界应力
/// σ_cr = π²E / λ²
pub fn euler_critical_stress(e : Double, lambda : Double) -> Double {
  if lambda <= 0.0 {
    0.0
  } else {
    PI * PI * e / square(lambda)
  }
}

///|
/// 临界长细比 (欧拉公式适用的最小长细比)
/// λ_p = π * sqrt(E / σ_p)
/// sigma_p: 比例极限
pub fn critical_slenderness(e : Double, sigma_p : Double) -> Double {
  if sigma_p <= 0.0 {
    0.0
  } else {
    PI * (e / sigma_p).sqrt()
  }
}

///|
/// 折减系数法计算稳定承载力
/// N = φ * A * f
/// phi: 稳定系数, a: 截面积, f: 设计强度
pub fn stability_capacity(phi : Double, area : Double, f : Double) -> Double {
  phi * area * f
}

///|
/// 查表或计算稳定系数 (简化经验公式)
/// 对于 Q235 钢，当 λ <= 100 时
/// φ ≈ 1 - 0.00668λ
/// 当 100 < λ <= 150 时
/// φ ≈ 1.1 - 0.0107λ (近似)
pub fn stability_factor_steel_q235(lambda : Double) -> Double {
  if lambda <= 0.0 {
    1.0
  } else if lambda <= 100.0 {
    max(0.1, 1.0 - 0.00668 * lambda)
  } else if lambda <= 150.0 {
    max(0.1, 1.1 - 0.0107 * lambda)
  } else {
    0.1 // 超出范围，返回最小值
  }
}

// ==================== 刚架分析 ====================

///|
/// 门式刚架侧移刚度 (两端固定柱)
/// K = 24EI_c / h³
/// i_c: 柱惯性矩, h: 柱高
pub fn portal_frame_lateral_stiffness_fixed(
  e : Double,
  i_c : Double,
  h : Double,
) -> Double {
  if h <= 0.0 {
    0.0
  } else {
    24.0 * e * i_c / cube(h)
  }
}

///|
/// 门式刚架侧移刚度 (底部铰接柱)
/// K = 3EI_c / h³
pub fn portal_frame_lateral_stiffness_hinged(
  e : Double,
  i_c : Double,
  h : Double,
) -> Double {
  if h <= 0.0 {
    0.0
  } else {
    3.0 * e * i_c / cube(h)
  }
}

///|
/// 刚架柱的抗侧刚度贡献
/// D = 12EI / (μ²h³)
/// mu: 柱的计算长度系数
pub fn column_lateral_stiffness(
  e : Double,
  i : Double,
  h : Double,
  mu : Double,
) -> Double {
  if h <= 0.0 || mu <= 0.0 {
    0.0
  } else {
    12.0 * e * i / (mu * mu * cube(h))
  }
}

// ==================== 虚功法 ====================

///|
/// 单位荷载法求位移
/// Δ = Σ(M * M̄ * L) / (E * I)
/// 这里提供单段的贡献计算
pub fn virtual_work_bending(
  m_actual : Double,
  m_virtual : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  if e <= 0.0 || i <= 0.0 {
    0.0
  } else {
    m_actual * m_virtual * l / (e * i)
  }
}

///|
/// 单位荷载法求位移 - 轴力贡献
/// Δ_N = Σ(N * N̄ * L) / (E * A)
pub fn virtual_work_axial(
  n_actual : Double,
  n_virtual : Double,
  l : Double,
  e : Double,
  a : Double,
) -> Double {
  if e <= 0.0 || a <= 0.0 {
    0.0
  } else {
    n_actual * n_virtual * l / (e * a)
  }
}

///|
/// 单位荷载法求位移 - 剪力贡献
/// Δ_V = Σ(κ * V * V̄ * L) / (G * A)
pub fn virtual_work_shear(
  v_actual : Double,
  v_virtual : Double,
  l : Double,
  g : Double,
  a : Double,
  kappa : Double,
) -> Double {
  if g <= 0.0 || a <= 0.0 {
    0.0
  } else {
    kappa * v_actual * v_virtual * l / (g * a)
  }
}

// ==================== 卡氏定理 ====================

///|
/// 卡氏第二定理 - 弹性体系位移
/// Δᵢ = ∂U/∂Pᵢ
/// 这里计算单根杆件的应变能对某力的偏导贡献
/// 弯矩: ∂U/∂P = ∫(M/EI)(∂M/∂P)dx
pub fn castigliano_bending_contribution(
  m : Double,
  dm_dp : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  if e <= 0.0 || i <= 0.0 {
    0.0
  } else {
    m * dm_dp * l / (e * i)
  }
}

// ==================== 连续梁分析 ====================

///|
/// 三弯矩方程系数
/// 对于等截面连续梁: M_{n-1}*L_n + 2*M_n*(L_n + L_{n+1}) + M_{n+1}*L_{n+1} = -6*Σ(荷载项)
/// 返回系数 (a, b, c) 即 a*M_{n-1} + b*M_n + c*M_{n+1} = rhs
pub fn three_moment_coefficients(
  l_left : Double,
  l_right : Double,
) -> (Double, Double, Double) {
  (l_left, 2.0 * (l_left + l_right), l_right)
}

///|
/// 三弯矩方程右端项 - 均布荷载
/// 左跨均布荷载 q_l，右跨均布荷载 q_r
/// rhs = -6 * (q_l*L_l³/24/L_l + q_r*L_r³/24/L_r) = -(q_l*L_l² + q_r*L_r²)/4
pub fn three_moment_rhs_uniform(
  q_left : Double,
  l_left : Double,
  q_right : Double,
  l_right : Double,
) -> Double {
  -(q_left * square(l_left) + q_right * square(l_right)) / 4.0
}

///|
/// 两跨连续梁支座弯矩 (两端简支，等跨等载)
/// 中间支座弯矩 MB = -qL²/8
pub fn continuous_beam_2span_uniform_mid_moment(
  q : Double,
  l : Double,
) -> Double {
  -q * square(l) / 8.0
}

///|
/// 三跨连续梁支座弯矩 (两端简支，等跨等载)
/// 简化公式: MB = MC = -qL²/10
pub fn continuous_beam_3span_uniform_support_moment(
  q : Double,
  l : Double,
) -> Double {
  -q * square(l) / 10.0
}

// ==================== 拱结构 ====================

///|
/// 三铰拱的水平推力 (竖向均布荷载)
/// H = qL² / (8f)
/// f: 拱的矢高
pub fn three_hinged_arch_thrust_uniform(
  q : Double,
  l : Double,
  f : Double,
) -> Double {
  if f <= 0.0 {
    0.0
  } else {
    q * square(l) / (8.0 * f)
  }
}

///|
/// 三铰拱的水平推力 (跨中集中力)
/// H = PL / (4f)
pub fn three_hinged_arch_thrust_point_mid(
  p : Double,
  l : Double,
  f : Double,
) -> Double {
  if f <= 0.0 {
    0.0
  } else {
    p * l / (4.0 * f)
  }
}

///|
/// 抛物线拱轴线方程
/// y = 4fx(L-x) / L²
pub fn parabolic_arch_y(x : Double, l : Double, f : Double) -> Double {
  if l <= 0.0 {
    0.0
  } else {
    4.0 * f * x * (l - x) / square(l)
  }
}

///|
/// 抛物线拱轴线斜率
/// dy/dx = 4f(L - 2x) / L²
pub fn parabolic_arch_slope(x : Double, l : Double, f : Double) -> Double {
  if l <= 0.0 {
    0.0
  } else {
    4.0 * f * (l - 2.0 * x) / square(l)
  }
}

///|
/// 圆弧拱的半径
/// R = (L²/4 + f²) / (2f)
pub fn circular_arch_radius(l : Double, f : Double) -> Double {
  if f <= 0.0 {
    0.0
  } else {
    (square(l) / 4.0 + square(f)) / (2.0 * f)
  }
}

// ==================== 应变能 ====================

///|
/// 弯曲应变能 (均匀弯矩段)
/// U = M²L / (2EI)
pub fn bending_strain_energy(
  m : Double,
  l : Double,
  e : Double,
  i : Double,
) -> Double {
  if e <= 0.0 || i <= 0.0 {
    0.0
  } else {
    square(m) * l / (2.0 * e * i)
  }
}

///|
/// 轴向应变能
/// U = N²L / (2EA)
pub fn axial_strain_energy(
  n : Double,
  l : Double,
  e : Double,
  a : Double,
) -> Double {
  if e <= 0.0 || a <= 0.0 {
    0.0
  } else {
    square(n) * l / (2.0 * e * a)
  }
}

///|
/// 剪切应变能
/// U = κV²L / (2GA)
pub fn shear_strain_energy(
  v : Double,
  l : Double,
  g : Double,
  a : Double,
  kappa : Double,
) -> Double {
  if g <= 0.0 || a <= 0.0 {
    0.0
  } else {
    kappa * square(v) * l / (2.0 * g * a)
  }
}

///|
/// 扭转应变能
/// U = T²L / (2GIp)
pub fn torsion_strain_energy(
  t : Double,
  l : Double,
  g : Double,
  ip : Double,
) -> Double {
  if g <= 0.0 || ip <= 0.0 {
    0.0
  } else {
    square(t) * l / (2.0 * g * ip)
  }
}

// ==================== 动力学基础 ====================

///|
/// 单自由度系统自振圆频率
/// ω = sqrt(k/m)
pub fn natural_frequency(k : Double, m : Double) -> Double {
  if m <= 0.0 || k <= 0.0 {
    0.0
  } else {
    (k / m).sqrt()
  }
}

///|
/// 单自由度系统自振周期
/// T = 2π/ω = 2π*sqrt(m/k)
pub fn natural_period(k : Double, m : Double) -> Double {
  if m <= 0.0 || k <= 0.0 {
    0.0
  } else {
    2.0 * PI * (m / k).sqrt()
  }
}

///|
/// 简支梁一阶自振频率 (集中质量在跨中)
/// ω₁ = sqrt(48EI / (mL³))
pub fn beam_frequency_concentrated_mass(
  e : Double,
  i : Double,
  m : Double,
  l : Double,
) -> Double {
  if m <= 0.0 || l <= 0.0 || e <= 0.0 || i <= 0.0 {
    0.0
  } else {
    (48.0 * e * i / (m * cube(l))).sqrt()
  }
}

///|
/// 简支梁一阶自振频率 (均布质量)
/// ω₁ = π²*sqrt(EI / (ρAL⁴))
/// 其中 m_total = ρ*A*L
pub fn beam_frequency_distributed_mass(
  e : Double,
  i : Double,
  rho_a : Double,
  l : Double,
) -> Double {
  if rho_a <= 0.0 || l <= 0.0 || e <= 0.0 || i <= 0.0 {
    0.0
  } else {
    PI * PI * (e * i / (rho_a * pow4(l))).sqrt()
  }
}

///|
/// 悬臂梁一阶自振频率 (端部集中质量)
/// ω₁ = sqrt(3EI / (mL³))
pub fn cantilever_frequency_tip_mass(
  e : Double,
  i : Double,
  m : Double,
  l : Double,
) -> Double {
  if m <= 0.0 || l <= 0.0 || e <= 0.0 || i <= 0.0 {
    0.0
  } else {
    (3.0 * e * i / (m * cube(l))).sqrt()
  }
}

///|
/// 阻尼比
/// ξ = c / (2*sqrt(k*m)) = c / (2*m*ω)
pub fn damping_ratio(c : Double, k : Double, m : Double) -> Double {
  if m <= 0.0 || k <= 0.0 {
    0.0
  } else {
    c / (2.0 * (k * m).sqrt())
  }
}

///|
/// 有阻尼自振频率
/// ω_d = ω * sqrt(1 - ξ²)
pub fn damped_frequency(omega : Double, xi : Double) -> Double {
  if xi >= 1.0 {
    0.0 // 过阻尼
  } else {
    omega * (1.0 - square(xi)).sqrt()
  }
}

// ==================== 组合结构 ====================

///|
/// 换算截面法 - 钢筋混凝土组合截面
/// 换算面积 A₀ = Ac + αₑ*As
/// alpha_e = Es/Ec
pub fn transformed_area(ac : Double, as_ : Double, alpha_e : Double) -> Double {
  ac + alpha_e * as_
}

///|
/// 换算惯性矩 (近似，钢筋在受拉侧)
/// I₀ = Ic + αₑ*As*(d-x)²
pub fn transformed_inertia(
  ic : Double,
  as_ : Double,
  alpha_e : Double,
  d : Double,
  x : Double,
) -> Double {
  ic + alpha_e * as_ * square(d - x)
}
