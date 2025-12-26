///|
fn main {
  // 示例：简支梁在均布荷载下的最大弯矩
  // q = 10 kN/m, L = 6 m
  let m_max = @lib.simply_supported_max_moment_uniform(10.0, 6.0)
  println("简支梁最大弯矩: \{m_max} kN·m")
  
  // 示例：欧拉临界力
  // E = 200 GPa, I = 0.0001 m⁴, L = 5 m, μ = 1
  let pcr = @lib.euler_critical_force(200.0e9, 0.0001, 5.0, 1.0)
  println("欧拉临界力: \{pcr / 1.0e6} MN")
}
