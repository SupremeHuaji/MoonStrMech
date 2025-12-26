// MoonStrMech - 工具函数

// ==================== 基本数学运算 ====================

///|
/// 计算平方
pub fn square(x : Double) -> Double {
  x * x
}

///|
/// 计算立方
pub fn cube(x : Double) -> Double {
  x * x * x
}

///|
/// 计算四次方
pub fn pow4(x : Double) -> Double {
  let x2 = x * x
  x2 * x2
}

///|
/// 计算绝对值
pub fn abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
/// 计算最大值
pub fn max(a : Double, b : Double) -> Double {
  if a > b { a } else { b }
}

///|
/// 计算最小值
pub fn min(a : Double, b : Double) -> Double {
  if a < b { a } else { b }
}

///|
/// 安全除法，分母为零时返回默认值
pub fn safe_div(num : Double, den : Double, fallback : Double) -> Double {
  if abs(den) <= EPSILON {
    fallback
  } else {
    num / den
  }
}

///|
/// 判断两个浮点数是否近似相等
pub fn approx_equal(a : Double, b : Double, tolerance : Double) -> Bool {
  abs(a - b) < tolerance
}

///|
/// 判断浮点数是否接近零
pub fn is_zero(x : Double) -> Bool {
  abs(x) < EPSILON
}

///|
/// 符号函数
pub fn sign(x : Double) -> Double {
  if x > EPSILON {
    1.0
  } else if x < -EPSILON {
    -1.0
  } else {
    0.0
  }
}

///|
/// 限制值在范围内
pub fn clamp(value : Double, min_val : Double, max_val : Double) -> Double {
  if value < min_val {
    min_val
  } else if value > max_val {
    max_val
  } else {
    value
  }
}

///|
/// 线性插值
pub fn lerp(a : Double, b : Double, t : Double) -> Double {
  a + (b - a) * t
}

///|
/// 角度转弧度
pub fn deg_to_rad(deg : Double) -> Double {
  deg * PI / 180.0
}

///|
/// 弧度转角度
pub fn rad_to_deg(rad : Double) -> Double {
  rad * 180.0 / PI
}

// ==================== 向量运算 (2D) ====================

///|
/// 二维向量点积
pub fn dot2d(x1 : Double, y1 : Double, x2 : Double, y2 : Double) -> Double {
  x1 * x2 + y1 * y2
}

///|
/// 二维向量叉积 (返回标量，z分量)
pub fn cross2d(x1 : Double, y1 : Double, x2 : Double, y2 : Double) -> Double {
  x1 * y2 - y1 * x2
}

///|
/// 二维向量长度
pub fn vec2d_length(x : Double, y : Double) -> Double {
  (x * x + y * y).sqrt()
}

///|
/// 两点间距离
pub fn distance2d(
  x1 : Double,
  y1 : Double,
  x2 : Double,
  y2 : Double
) -> Double {
  let dx = x2 - x1
  let dy = y2 - y1
  (dx * dx + dy * dy).sqrt()
}

// ==================== 简单矩阵运算 (2x2, 3x3) ====================

///|
/// 2x2矩阵行列式
/// | a b |
/// | c d |
pub fn det2x2(a : Double, b : Double, c : Double, d : Double) -> Double {
  a * d - b * c
}

///|
/// 3x3矩阵行列式
/// | a b c |
/// | d e f |
/// | g h i |
pub fn det3x3(
  a : Double,
  b : Double,
  c : Double,
  d : Double,
  e : Double,
  f : Double,
  g : Double,
  h : Double,
  i : Double
) -> Double {
  a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

///|
/// 2x2矩阵求逆，返回 (a', b', c', d') 或全零表示奇异
pub fn inv2x2(
  a : Double,
  b : Double,
  c : Double,
  d : Double
) -> (Double, Double, Double, Double) {
  let det = det2x2(a, b, c, d)
  if abs(det) < EPSILON {
    (0.0, 0.0, 0.0, 0.0)
  } else {
    let inv_det = 1.0 / det
    (d * inv_det, -b * inv_det, -c * inv_det, a * inv_det)
  }
}

///|
/// 2x2矩阵乘向量 [a b; c d] * [x; y]
pub fn mat2x2_vec(
  a : Double,
  b : Double,
  c : Double,
  d : Double,
  x : Double,
  y : Double
) -> (Double, Double) {
  (a * x + b * y, c * x + d * y)
}

///|
/// 求解2x2线性方程组 [a b; c d] * [x; y] = [e; f]
pub fn solve2x2(
  a : Double,
  b : Double,
  c : Double,
  d : Double,
  e : Double,
  f : Double
) -> (Double, Double) {
  let det = det2x2(a, b, c, d)
  if abs(det) < EPSILON {
    (0.0, 0.0) // 奇异矩阵
  } else {
    let x = (e * d - b * f) / det
    let y = (a * f - e * c) / det
    (x, y)
  }
}

///|
/// 求解3x3线性方程组 (克莱姆法则)
/// 返回 (x, y, z)
pub fn solve3x3(
  a11 : Double,
  a12 : Double,
  a13 : Double,
  a21 : Double,
  a22 : Double,
  a23 : Double,
  a31 : Double,
  a32 : Double,
  a33 : Double,
  b1 : Double,
  b2 : Double,
  b3 : Double
) -> (Double, Double, Double) {
  let det = det3x3(a11, a12, a13, a21, a22, a23, a31, a32, a33)
  if abs(det) < EPSILON {
    (0.0, 0.0, 0.0)
  } else {
    let det_x = det3x3(b1, a12, a13, b2, a22, a23, b3, a32, a33)
    let det_y = det3x3(a11, b1, a13, a21, b2, a23, a31, b3, a33)
    let det_z = det3x3(a11, a12, b1, a21, a22, b2, a31, a32, b3)
    (det_x / det, det_y / det, det_z / det)
  }
}

// ==================== 数值积分 ====================

///|
/// 梯形积分法 (单段)
/// 积分 f 在 [a, b] 上的值，给定 f(a) 和 f(b)
pub fn trapezoid_integral(fa : Double, fb : Double, length : Double) -> Double {
  (fa + fb) * length / 2.0
}

///|
/// 辛普森积分 (单段)
/// 积分 f 在 [a, b] 上的值，给定 f(a), f(mid), f(b)
pub fn simpson_integral(
  fa : Double,
  fmid : Double,
  fb : Double,
  length : Double
) -> Double {
  (fa + 4.0 * fmid + fb) * length / 6.0
}

// ==================== 多项式运算 ====================

///|
/// 一次多项式求值 y = a*x + b
pub fn poly1_eval(a : Double, b : Double, x : Double) -> Double {
  a * x + b
}

///|
/// 二次多项式求值 y = a*x² + b*x + c
pub fn poly2_eval(a : Double, b : Double, c : Double, x : Double) -> Double {
  (a * x + b) * x + c
}

///|
/// 三次多项式求值 y = a*x³ + b*x² + c*x + d
pub fn poly3_eval(
  a : Double,
  b : Double,
  c : Double,
  d : Double,
  x : Double
) -> Double {
  ((a * x + b) * x + c) * x + d
}

// ==================== 几何工具 ====================

///|
/// 计算三角形面积 (给定三边长，海伦公式)
pub fn triangle_area_heron(a : Double, b : Double, c : Double) -> Double {
  let s = (a + b + c) / 2.0
  let area_sq = s * (s - a) * (s - b) * (s - c)
  if area_sq <= 0.0 {
    0.0
  } else {
    area_sq.sqrt()
  }
}

///|
/// 计算三角形面积 (给定三点坐标)
pub fn triangle_area_coords(
  x1 : Double,
  y1 : Double,
  x2 : Double,
  y2 : Double,
  x3 : Double,
  y3 : Double
) -> Double {
  abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0)
}

///|
/// 计算多边形面积 (顶点数组，顺或逆时针，Shoelace公式)
pub fn polygon_area(coords : Array[(Double, Double)]) -> Double {
  let n = coords.length()
  if n < 3 {
    return 0.0
  }
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let (x1, y1) = coords[i]
    let next = if i + 1 < n { i + 1 } else { 0 }
    let (x2, y2) = coords[next]
    sum = sum + (x1 * y2 - x2 * y1)
  }
  abs(sum) / 2.0
}

///|
/// 计算多边形形心 (顶点数组)
pub fn polygon_centroid(
  coords : Array[(Double, Double)]
) -> (Double, Double) {
  let n = coords.length()
  if n < 3 {
    return (0.0, 0.0)
  }
  let mut cx = 0.0
  let mut cy = 0.0
  let mut area_sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let (x1, y1) = coords[i]
    let next = if i + 1 < n { i + 1 } else { 0 }
    let (x2, y2) = coords[next]
    let cross = x1 * y2 - x2 * y1
    area_sum = area_sum + cross
    cx = cx + (x1 + x2) * cross
    cy = cy + (y1 + y2) * cross
  }
  if abs(area_sum) < EPSILON {
    (0.0, 0.0)
  } else {
    let factor = 1.0 / (3.0 * area_sum)
    (cx * factor, cy * factor)
  }
}

// ==================== 坐标变换 ====================

///|
/// 二维旋转 (绕原点)
/// 返回旋转后的坐标 (x', y')
/// 使用泰勒级数近似 sin 和 cos
pub fn rotate2d(
  x : Double,
  y : Double,
  angle_rad : Double
) -> (Double, Double) {
  let cos_a = cos_approx(angle_rad)
  let sin_a = sin_approx(angle_rad)
  (x * cos_a - y * sin_a, x * sin_a + y * cos_a)
}

///|
/// 正弦函数近似 (泰勒级数)
pub fn sin_approx(x : Double) -> Double {
  // 归一化到 [-π, π]
  let mut a = x
  while a > PI {
    a = a - 2.0 * PI
  }
  while a < -PI {
    a = a + 2.0 * PI
  }
  // 泰勒级数: sin(x) ≈ x - x³/6 + x⁵/120 - x⁷/5040
  let x2 = a * a
  let x3 = x2 * a
  let x5 = x3 * x2
  let x7 = x5 * x2
  a - x3 / 6.0 + x5 / 120.0 - x7 / 5040.0
}

///|
/// 余弦函数近似 (泰勒级数)
pub fn cos_approx(x : Double) -> Double {
  // 归一化到 [-π, π]
  let mut a = x
  while a > PI {
    a = a - 2.0 * PI
  }
  while a < -PI {
    a = a + 2.0 * PI
  }
  // 泰勒级数: cos(x) ≈ 1 - x²/2 + x⁴/24 - x⁶/720
  let x2 = a * a
  let x4 = x2 * x2
  let x6 = x4 * x2
  1.0 - x2 / 2.0 + x4 / 24.0 - x6 / 720.0
}

///|
/// 二维平移
pub fn translate2d(
  x : Double,
  y : Double,
  dx : Double,
  dy : Double
) -> (Double, Double) {
  (x + dx, y + dy)
}

// ==================== 统计工具 ====================

///|
/// 数组求和
pub fn array_sum(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  for i = 0; i < arr.length(); i = i + 1 {
    sum = sum + arr[i]
  }
  sum
}

///|
/// 数组平均值
pub fn array_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    0.0
  } else {
    array_sum(arr) / arr.length().to_double()
  }
}

///|
/// 数组最大值
pub fn array_max(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut max_val = arr[0]
  for i = 1; i < arr.length(); i = i + 1 {
    if arr[i] > max_val {
      max_val = arr[i]
    }
  }
  max_val
}

///|
/// 数组最小值
pub fn array_min(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut min_val = arr[0]
  for i = 1; i < arr.length(); i = i + 1 {
    if arr[i] < min_val {
      min_val = arr[i]
    }
  }
  min_val
}

