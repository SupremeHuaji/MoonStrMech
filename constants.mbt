// MoonStrMech - 结构力学常量定义

// ==================== 数学常量 ====================

///|
/// 圆周率
pub const PI : Double = 3.141592653589793

///|
/// 自然对数底数
pub const E_CONST : Double = 2.718281828459045

///|
/// 浮点数比较精度
pub const EPSILON : Double = 1.0e-10

// ==================== 长度系数 (压杆稳定) ====================
// 不同端部约束条件下的长度系数 μ

///|
/// 两端铰支 (Hinged-Hinged)
pub const MU_HH : Double = 1.0

///|
/// 一端固定一端自由 (Fixed-Free) 悬臂
pub const MU_FF : Double = 2.0

///|
/// 两端固定 (Fixed-Fixed)
pub const MU_FIXED : Double = 0.5

///|
/// 一端固定一端铰支 (Fixed-Hinged)
pub const MU_FH : Double = 0.7

// ==================== 常用材料弹性模量 (GPa) ====================

///|
/// 钢材弹性模量
pub const E_STEEL : Double = 206.0

///|
/// 混凝土弹性模量 (C30)
pub const E_CONCRETE_C30 : Double = 30.0

///|
/// 混凝土弹性模量 (C40)
pub const E_CONCRETE_C40 : Double = 32.5

///|
/// 木材弹性模量
pub const E_WOOD : Double = 10.0

///|
/// 铝合金弹性模量
pub const E_ALUMINUM : Double = 70.0

// ==================== 常用材料泊松比 ====================

///|
pub const NU_STEEL : Double = 0.3

///|
pub const NU_CONCRETE : Double = 0.2

///|
pub const NU_WOOD : Double = 0.35

///|
pub const NU_ALUMINUM : Double = 0.33

// ==================== 常用材料密度 (kg/m³) ====================

///|
pub const RHO_STEEL : Double = 7850.0

///|
pub const RHO_CONCRETE : Double = 2500.0

///|
pub const RHO_WOOD : Double = 600.0

///|
pub const RHO_ALUMINUM : Double = 2700.0

// ==================== 常用荷载 (kN/m³) ====================

///|
/// 钢筋混凝土容重
pub const GAMMA_RC : Double = 25.0

///|
/// 钢材容重
pub const GAMMA_STEEL : Double = 78.5

///|
/// 木材容重
pub const GAMMA_WOOD : Double = 6.0

///|
/// 水容重
pub const GAMMA_WATER : Double = 10.0

///|
/// 土容重 (一般粘性土)
pub const GAMMA_SOIL : Double = 18.0

// ==================== 截面形状系数 (剪切) ====================

///|
/// 矩形截面剪切形状系数
pub const KAPPA_RECT : Double = 1.2

///|
/// 圆形截面剪切形状系数
pub const KAPPA_CIRCLE : Double = 1.11

///|
/// 工字形截面剪切形状系数 (近似)
pub const KAPPA_I_SECTION : Double = 1.0

// ==================== 力矩分配法相关 ====================

///|
/// 远端铰支传递系数
pub const CARRYOVER_HINGED : Double = 0.0

///|
/// 远端固定传递系数
pub const CARRYOVER_FIXED : Double = 0.5

// ==================== 单位换算 ====================

///|
/// kN 转 N
pub const KN_TO_N : Double = 1000.0

///|
/// MPa 转 Pa
pub const MPA_TO_PA : Double = 1000000.0

///|
/// GPa 转 Pa
pub const GPA_TO_PA : Double = 1000000000.0

///|
/// m 转 mm
pub const M_TO_MM : Double = 1000.0

///|
/// mm 转 m
pub const MM_TO_M : Double = 0.001

