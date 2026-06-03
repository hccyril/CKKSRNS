"""
FullRNS-HEAAN / CKKS 教学 demo
================================

本文件是论文
    "A Full RNS Variant of Approximate Homomorphic Encryption"
    (https://eprint.iacr.org/2018/931, SAC 2018)
的 Python 教学演示实现。

目标：用尽量小、尽量可读的 Python 代码，将论文中的核心方案步骤完整串起来。
**不追求安全性或性能**，仅用于学习。

===== 本仓库 C++ 代码与论文的对应关系 =====

    C++ 模块              论文章节           功能
    ────────────────────────────────────────────────────────────
    src/Context.*         论文 2.1-2.3       环运算、NTT、RNS 基转换、raise/back/rescale
    src/Scheme.*          论文 4             KeyGen、加解密、同态加/减/乘、重缩放、旋转、共轭
    src/SecretKey.*       论文 4 KeyGen      小汉明重量(HWT) secret key 采样
    src/Ciphertext.*      论文 4             密文数据结构
    src/Key.*             论文 4 KSGen       密钥开关(key switching)密钥
    src/SchemeAlgo.*      论文 5             高阶算法(多项式求值等，本 demo 不涉及)

===== 本 demo 与论文各节的对应关系 =====

    论文 2.1  : CKKS/HEAAN 近似数编码/解码、缩放因子(scale)概念、同态加乘与 rescale 思想
    论文 2.2  : RNS 表示 [a]_B —— 多项式系数按多个模数分别存储
    论文 2.3  : Fast Basis Conversion Conv_{C->B} / Conv_{B->C}
    论文 3.1  : approximate basis —— 各 q_j 近似同一个固定 base q
    论文 3.2  : Algorithm 1 ModUp_{C->D}, Algorithm 2 ModDown_{D->C}
    论文 4    : Setup, KSGen, KeyGen, Enc, Dec, Add, Sub, Mult, RS,
                以及 conjugation key / rotation key 的生成与使用

===== 与 C++ 实现的差异说明 =====

    1. 多项式表示：
       真实 C++ 库把每个 RNS 分量存为 NTT(数论变换)表示，做快速逐点乘法；
       本 demo 把多项式保存在"系数表示"中，并逐模数做 O(N^2) 环乘法。
       这不改变 RNS 层面的公式，只是去掉了 NTT 优化层，便于阅读。

    2. 乘法 Karatsuba 技巧：
       C++ Scheme::mult 用 Karatsuba 恒等式
           u1 = (c0+c1)(d0+d1) - c0*d0 - c1*d1
       代替直接计算 c0*d1 + c1*d0，减少一次环乘。
       本 demo 直接使用三乘积公式(u0, u1, u2)，更直观。
       两者数学等价。

    3. Eval key 中的 P 因子：
       C++ addMultKey 对 s^2 先做 evalAndEqual (即乘 P) 再混入 evk；
       本 demo 在 evk 构造中直接写 P*s^2，效果等价。

    4. 参数极小：N=8，slots=N/2=4 个复数。
       这完全没有安全性，误差也会明显偏大——尤其是乘法+重线性化+rescale 之后。

    5. 随机性使用固定 seed，便于复现。
"""

from __future__ import annotations

from dataclasses import dataclass
import cmath
import math
import random
from typing import Dict, Iterable, List, Sequence, Tuple


# ============================================================================
# 类型别名
# ============================================================================
# Poly: 一个 degree < N 的多项式，系数用 Python int 存储（可能为负或超大）
Poly = List[int]
# RNSPoly: 论文 2.2 的 [a]_B 表示——同一个多项式对多个模数分别取余
#   key = 模数 q_j, value = 系数列表(每个系数 in [0, q_j))
RNSPoly = Dict[int, Poly]
# RNSTuple: 密文的两个分量 (c0, c1)
RNSTuple = Tuple[RNSPoly, RNSPoly]


# ============================================================================
# 1. 基础工具函数
# ============================================================================


def inv_mod(a: int, q: int) -> int:
    """
    返回 a^{-1} mod q。
    论文全篇频繁使用模逆运算，例如 basis conversion 中的 qhat_j^{-1} mod q_j。
    Python 3.8+ 内置 pow(a, -1, q)。
    """
    return pow(a % q, -1, q)


def centered(x: int, q: int) -> int:
    """
    论文第 2 节记号 [a]_q 的中心代表元：取 (-q/2, q/2] 中的整数。

    在 CKKS 解密时，我们需要把模 q 下的值还原为"小的"有符号整数，
    才能正确解读为带缩放因子的明文。
    """
    x %= q
    return x if x <= q // 2 else x - q


# ---------- 多项式环 R = Z[X]/(X^N + 1) 的基本运算 ----------
# 论文第 2 节: 所有运算都在环 R_q = Z_q[X]/(X^N+1) 中进行


def poly_zero(N: int) -> Poly:
    """返回零多项式。"""
    return [0] * N


def poly_add(a: Poly, b: Poly, q: int) -> Poly:
    """
    环 R_q 中的多项式加法: (a + b) mod q。
    论文第 4 节 Add/Sub 操作的基础。
    """
    return [(x + y) % q for x, y in zip(a, b)]


def poly_sub(a: Poly, b: Poly, q: int) -> Poly:
    """环 R_q 中的多项式减法: (a - b) mod q。"""
    return [(x - y) % q for x, y in zip(a, b)]


def poly_neg(a: Poly, q: int) -> Poly:
    """环 R_q 中的多项式取负: -a mod q。"""
    return [(-x) % q for x in a]


def poly_scalar_mul(a: Poly, s: int, q: int) -> Poly:
    """环 R_q 中的标量乘: s * a mod q。"""
    return [(x * s) % q for x in a]


def poly_mul(a: Poly, b: Poly, q: int) -> Poly:
    """
    环 R_q = Z_q[X]/(X^N + 1) 中的乘法（学校本卷积）。

    对普通卷积中次数 >= N 的项，使用 X^N ≡ -1 (mod X^N+1) 折回：
        X^{N+t} = -X^t

    论文第 4 节所有密文分量都属于某个 R_{q_j} 或 R_{p_i}。

    注意: C++ 代码在 NTT 表示下做逐点乘法，时间复杂度 O(N log N)；
    这里用 O(N^2) 的学校本卷积，便于理解。
    """
    N = len(a)
    out = [0] * N
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            deg = i + j
            if deg < N:
                out[deg] += ai * bj
            else:
                # X^{deg} = X^{N + (deg-N)} = -X^{deg-N}
                out[deg - N] -= ai * bj
    return [x % q for x in out]


def poly_mod(poly: Poly, q: int) -> Poly:
    """将多项式的每个系数模 q。"""
    return [x % q for x in poly]


def poly_conjugate(a: Poly) -> Poly:
    """
    论文共轭操作 (参见 C++ Context::conjugate)：

    在环 Z[X]/(X^N+1) 中，共轭映射为:
        a(X) -> a(X^{-1}) = a(X^{2N-1})

    因为 X^N = -1，所以:
        X^{-1} = X^{N-1} * X^{-N} = -X^{N-1}

    对于 a(X) = a_0 + a_1*X + ... + a_{N-1}*X^{N-1}：
        a(X^{-1}) = a_0 - a_{N-1}*X - a_{N-2}*X^2 - ... - a_1*X^{N-1}

    系数表示: [a_0, a_1, ..., a_{N-1}] -> [a_0, -a_{N-1}, -a_{N-2}, ..., -a_1]

    这等价于 C++ 中的:
        resi[n] = ai[N - 1 - n]  (配合 NTT 域的特殊存储顺序)
    """
    N = len(a)
    result = [a[0]]
    for i in range(1, N):
        result.append(-a[N - i])
    return result


def poly_rotate(a: Poly, power: int) -> Poly:
    """
    论文左旋转操作 (参见 C++ Context::leftRot)：

    在环 Z[X]/(X^N+1) 中，将 a(X) 映射为 a(X^power)。

    对于 a(X) = sum_n a_n * X^n：
        a(X^power) = sum_n a_n * X^{n*power}

    由于 X^N = -1，当 n*power mod 2N >= N 时：
        X^{n*power mod 2N} = -X^{(n*power mod 2N) - N}

    参数:
        a:     多项式系数
        power: 旋转幂次，由 rotGroup 决定 (如 5^r mod 2N)
    """
    N = len(a)
    M = 2 * N
    result = [0] * N
    for n in range(N):
        shift = (n * power) % M
        if shift < N:
            result[shift] += a[n]
        else:
            # X^{shift} = X^{N + (shift-N)} = -X^{shift-N}
            result[shift - N] -= a[n]
    return result


# ---------- RNS 层面的运算 ----------
# 论文 2.2: RNS 表示 [a]_B = (a mod q_0, a mod q_1, ..., a mod q_l)


def rns_from_poly(poly: Poly, moduli: Sequence[int]) -> RNSPoly:
    """
    论文 2.2 的 RNS 表示 [a]_B:
    将整数多项式按每个模数 q_j 取余，得到各分量。
    """
    return {q: poly_mod(poly, q) for q in moduli}


def rns_add(a: RNSPoly, b: RNSPoly) -> RNSPoly:
    """RNS 逐分量加法。论文第 4 节 Add 的基础。"""
    return {q: poly_add(a[q], b[q], q) for q in a}


def rns_sub(a: RNSPoly, b: RNSPoly) -> RNSPoly:
    """RNS 逐分量减法。论文第 4 节 Sub 的基础。"""
    return {q: poly_sub(a[q], b[q], q) for q in a}


def rns_neg(a: RNSPoly) -> RNSPoly:
    """RNS 逐分量取负。"""
    return {q: poly_neg(a[q], q) for q in a}


def rns_mul(a: RNSPoly, b: RNSPoly) -> RNSPoly:
    """
    RNS 逐分量环乘法。
    论文第 4 节的 Mult、KSGen 等大量使用此操作。
    每个模数 q_j 上独立做 R_{q_j} 中的环乘。
    """
    return {q: poly_mul(a[q], b[q], q) for q in a}


def rns_scalar_mul(a: RNSPoly, s: int) -> RNSPoly:
    """RNS 标量乘: 每个分量乘同一个整数 s。"""
    return {q: poly_scalar_mul(a[q], s, q) for q in a}


def rns_restrict(a: RNSPoly, moduli: Sequence[int]) -> RNSPoly:
    """
    RNS 子集提取: 从 [a]_B 中取出 [a]_{B'}，B' ⊂ B。
    论文中密文在不同 level 使用不同子集 C_l = {q_0, ..., q_l}。
    """
    return {q: a[q][:] for q in moduli}


# ============================================================================
# 2. 论文 2.3: Fast Basis Conversion (快速基转换)
# ============================================================================


def product(xs: Iterable[int]) -> int:
    """计算可迭代对象的乘积。"""
    r = 1
    for x in xs:
        r *= x
    return r


def fast_convert(poly_in_basis: RNSPoly, src: Sequence[int], dst: Sequence[int]) -> RNSPoly:
    """
    论文 2.3 Fast Basis Conversion:

        Conv_{src -> dst}([a]_src)_i
          = sum_j ( [a^{(j)} * qhat_j^{-1}]_{q_j} * qhat_j ) mod p_i

    其中 qhat_j = prod_{t != j} q_t  (论文中的 hat{q}_j)。

    重要性质（论文 Lemma 1 / Section 2.3）:
        输出通常不是精确整数 a 在新基上的 residue，而是
        a + Q * e 的 residue，其中 e 较小。
        论文 3.2 的 ModUp/ModDown 正是利用这个"近似但足够好"的性质。

    C++ 对应: Context::raiseAndEqual / Context::backAndEqual 内部的核心计算。
    """
    N = len(next(iter(poly_in_basis.values())))
    src_prod = product(src)
    out = {p: [0] * N for p in dst}

    for j, qj in enumerate(src):
        # qhat_j = prod_{t != j} src_t
        qhat = src_prod // qj
        # [qhat_j^{-1}]_{q_j}
        qhat_inv_mod_qj = inv_mod(qhat, qj)
        # [a^{(j)} * qhat_j^{-1}]_{q_j} —— 每个系数乘以 qhat_j^{-1} mod q_j
        scaled_coeffs = [(c * qhat_inv_mod_qj) % qj for c in poly_in_basis[qj]]

        # 对每个目标模数 p_i，累加 scaled_coeffs * (qhat_j mod p_i)
        for p in dst:
            qhat_mod_p = qhat % p
            for n in range(N):
                out[p][n] = (out[p][n] + scaled_coeffs[n] * qhat_mod_p) % p
    return out


# ============================================================================
# 3. 论文 3.2: ModUp 与 ModDown (近似模数提升/降低)
# ============================================================================


def mod_up_C_to_D(a_C: RNSPoly, C: Sequence[int], B: Sequence[int]) -> RNSPoly:
    """
    论文 3.2 Algorithm 1: Approximate Modulus Raising ModUp_{C->D}

    输入:  [a]_C, 其中 C = {q_0, ..., q_l}
    输出:  [a']_D, 其中 D = B ∪ C, 且 a' ≈ a (mod Q_l), |a'| <~ P*Q_l

    步骤:
        1. 用 Fast Basis Conversion Conv_{C->B} 生成特殊基 B 上的 residue
        2. 拼接 B 分量和原 C 分量

    C++ 对应: Context::raiseAndEqual 中的 ModUp 部分。

    注意: a' 可能不等于 a，而是 a + Q_l * v 的形式（v 由 fast_convert 的
    近似性引入），但只要 v 足够小，后续运算仍可正确解密。
    """
    a_B = fast_convert(a_C, C, B)
    # 拼接: D = B ∪ C
    return {**a_B, **{q: a_C[q][:] for q in C}}


def mod_down_D_to_C(a_D: RNSPoly, C: Sequence[int], B: Sequence[int]) -> RNSPoly:
    """
    论文 3.2 Algorithm 2: Approximate Modulus Reduction ModDown_{D->C}

    目标:
        输入 [b]_D (D = B ∪ C)，输出 [b']_C，使 b' ≈ P^{-1} * b。
        这里 P = prod(B)。

    公式对应论文 Algorithm 2:
        1. a_C = Conv_{B->C}([b]_B)            // 从 B 分量恢复近似的 a
        2. b'_j = P^{-1} * (b_{C,j} - a_{C,j}) mod q_j   // 减 a 后除以 P

    直觉:
        B 分量编码了 b mod P 的信息。
        从 B 恢复一个"小的" a ≈ b (mod P)，
        于是 b - a 可被 P (近似)整除；
        在 C 分量上乘以 P^{-1} 就等价于除以 P。

    C++ 对应: Context::backAndEqual 的核心逻辑。
    """
    P = product(B)
    # 提取 B 分量
    a_B = {p: a_D[p] for p in B}
    # 从 B 基转换到 C 基，得到近似的 a 在 C 上的表示
    approx_a_C = fast_convert(a_B, B, C)
    out = {}
    for q in C:
        pinv = inv_mod(P, q)
        # b'_j = P^{-1} * (b_j - a_j) mod q_j
        out[q] = [((a_D[q][i] - approx_a_C[q][i]) * pinv) % q
                  for i in range(len(a_D[q]))]
    return out


# ============================================================================
# 4. 论文 2.1: CKKS 编码 / 解码 (canonical embedding 的逆/正)
# ============================================================================


def primitive_root_exponents_for_slots(N: int) -> List[int]:
    """
    论文 2.1 使用 canonical embedding 的一半坐标。

    对于分圆多项式 Φ_M(X) = X^N + 1 (M = 2N)，其根为 ζ^k (k 为奇数)，
    其中 ζ = exp(2πi / M) = exp(πi / N) 是 2N 次本原单位根。

    CKKS 选取一组共轭代表元:
        {ζ^1, ζ^5, ζ^9, ..., ζ^{2N-3}}  (即 ζ^{4j+1}, j=0,...,N/2-1)

    对 N=8: 指数为 [1, 5, 9, 13]，共 slots = N/2 = 4 个 slot。

    C++ 对应: Context 中的 rotGroup 和 ksiPows 预计算。
    """
    return [(1 + 4 * j) % (2 * N) for j in range(N // 2)]


def encode(values: Sequence[complex], N: int, scale: int) -> Poly:
    """
    论文 2.1 的 Encode = σ^{-1} 后乘 scale 并舍入到 R。

    给定 slots = N/2 个复数 z_0, ..., z_{slots-1}，
    要找多项式 m(X) (degree < N) 使得:
        m(ζ^{e_j}) ≈ scale * z_j,  j = 0, ..., slots-1

    其中 e_j 是 primitive_root_exponents_for_slots 返回的指数。

    为保证 m(X) 有实系数(近似)，对未显式存储的共轭坐标填入 conjugate(z_j):
        m(ζ^{-e_j}) = conjugate(m(ζ^{e_j})) = scale * conjugate(z_j)

    最后用逆 DFT 风格的公式恢复 degree < N 的系数:
        c_k = (1/N) * sum_{e in odd_exponents} y_e * ζ^{-e*k}

    C++ 对应: Context::encode 中的 fftSpecialInv + 缩放 + 系数映射。
    本 demo 用直接的矩阵求逆代替 special FFT，牺牲效率换可读性。

    注意: 由于 scale 有限且多项式系数必须是整数，编码本身就会引入舍入误差。
    """
    slots = N // 2
    if len(values) != slots:
        raise ValueError(f"N={N} 时需要 {slots} 个 slot 值，但收到 {len(values)} 个")

    # ζ = exp(πi / N) 是 2N 次本原单位根
    zeta = cmath.exp(1j * cmath.pi / N)
    selected = primitive_root_exponents_for_slots(N)

    # 构造所有 2N 次奇数单位根处的目标值
    # 对选中的 e: y(ζ^e) = scale * z_j
    # 对其共轭 2N-e: y(ζ^{2N-e}) = scale * conjugate(z_j)
    y_by_exp = {}
    for e, z in zip(selected, values):
        y_by_exp[e] = scale * z
        y_by_exp[(-e) % (2 * N)] = scale * z.conjugate()

    # 逆 DFT: c_k = (1/N) * sum_e y(ζ^e) * ζ^{-e*k}
    # 这里 odd_exponents 覆盖所有 N 个奇数 e in {1, 3, 5, ..., 2N-1}
    odd_exponents = [2 * r + 1 for r in range(N)]
    coeffs = []
    for k in range(N):
        ck = sum(y_by_exp[e] * (zeta ** (-e * k)) for e in odd_exponents) / N
        coeffs.append(int(round(ck.real)))
    return coeffs


def decode(poly: Poly, N: int, scale: int) -> List[complex]:
    """
    论文 2.1 的 Decode = σ(m) / scale。

    将多项式在选定的 N/2 个单位根处求值，再除以 scale 恢复复数。

    对于第 j 个 slot:
        z_j = m(ζ^{e_j}) / scale

    C++ 对应: Context::decode 中的 qiINTT + fftSpecial + 除以 p。
    """
    zeta = cmath.exp(1j * cmath.pi / N)
    out = []
    for e in primitive_root_exponents_for_slots(N):
        # m(ζ^e) = sum_{k=0}^{N-1} poly[k] * ζ^{e*k}
        val = sum(poly[k] * (zeta ** (e * k)) for k in range(N)) / scale
        out.append(val)
    return out


# ============================================================================
# 5. 参数、密钥、密文数据结构：论文第 4 节 Setup
# ============================================================================


@dataclass
class Params:
    """
    论文第 4 节 Setup(q, L, ε; λ) 的教学参数。

    论文要求:
        - N 是 2 的幂，环 R = Z[X]/(X^N+1)
        - q_j ≈ 固定 base q (论文 3.1 approximate basis)
        - q_j ≡ 1 (mod 2N) 以支持 NTT (本 demo 不使用 NTT，但仍保持此约束)
        - B = {p_0, ..., p_{K-1}} 是特殊模数集，用于 key switching
        - P = prod(B) 要足够大，使 ModDown 的近似误差相对较小

    本 demo:
        N = 8 (所以 2N = 16, slots = 4)
        L = 2 (最大密文 level，支持 2 层乘法后 rescale)
        q_i 和 p_i 都选择 ≡ 1 (mod 16) 的小素数
        scale = 2^8 = 256 (缩放因子，论文中记为 Δ)

    注意: 这些参数完全没有安全性！仅用于教学演示。
    """

    N: int = 8                        # 环维度: R = Z[X]/(X^8 + 1)
    L: int = 2                        # 最大 level (密文可用模数链长度 = L+1)
    scale: int = 2**8                 # 缩放因子 Δ = 256
    # 密文模数链 C = {q_0, q_1, q_2}，level l 的密文使用 C_l = {q_0, ..., q_l}
    q: Tuple[int, ...] = (40961, 65537, 114689)
    # 特殊模数集 B = {p_0, p_1}，用于 key switching 中的 ModUp/ModDown
    # P = p_0 * p_1 = 147457 * 163841 ≈ 2^34
    B: Tuple[int, ...] = (147457, 163841)

    def C(self, level: int) -> Tuple[int, ...]:
        """返回 level 对应的密文模数集 C_l = {q_0, ..., q_l}。"""
        return self.q[: level + 1]

    def D(self, level: int) -> Tuple[int, ...]:
        """返回扩展基 D_l = B ∪ C_l，用于 key switching 中间步骤。"""
        return self.B + self.C(level)


@dataclass
class SecretKey:
    """
    论文第 4 节 sk = (1, s)。
    s 是小多项式(系数来自 {-1, 0, 1})。
    C++ 中用 sampleHWT (汉明重量采样)。
    """
    s: Poly


@dataclass
class PublicKey:
    """
    论文第 4 节 pk = (b, a) ∈ R_Q^2。
    满足 b + a*s ≈ e (mod Q)，即 pk 是"近似零"的加密。
    """
    b: RNSPoly
    a: RNSPoly


@dataclass
class EvalKey:
    """
    论文第 4 节 KSGen 的输出 evk = (b, a) ∈ R_{P*Q}^2。
    满足 b + a*s ≈ P * s' + e (mod P*Q)，
    其中 s' 是 key switching 的目标秘密(如 s^2, s_conj, s_rot)。
    """
    b: RNSPoly
    a: RNSPoly


@dataclass
class Ciphertext:
    """
    论文第 4 节密文 ct = (c0, c1) ∈ R_{Q_l}^2。
    解密: c0 + c1 * s ≈ Δ * m (mod Q_l)。
    level 表示当前使用的模数链层数 (C_level)。
    """
    c0: RNSPoly
    c1: RNSPoly
    level: int
    slots: int


# ============================================================================
# 6. FullRNS-CKKS 方案核心实现：论文第 4 节
# ============================================================================


class FullRNSCKKSDemo:
    """
    论文 "A Full RNS Variant of Approximate Homomorphic Encryption" 的完整教学实现。

    包含:
        - Setup / KeyGen / KSGen (论文 4.1)
        - Enc / Dec (论文 4.2)
        - Add / Sub (论文 4.3)
        - Mult / RS (论文 4.4)
        - multByConst / addConst (论文 4 常数操作)
        - conjugate / leftRotate (论文 4 槽位操作)
    """

    def __init__(self, params: Params, seed: int = 7):
        """
        论文第 4 节 Setup: 初始化参数和伪随机数生成器。
        使用固定 seed 保证结果可复现。
        """
        self.par = params
        self.rng = random.Random(seed)

    # ======================== 采样分布 ========================
    # 论文第 4 节 Setup 定义 key, enc, err 分布
    # 附录 A 讨论噪声分析
    # 这里使用极小噪声(系数 ∈ {-1,0,1})，方便肉眼观察

    def sample_secret(self) -> Poly:
        """
        论文 KeyGen 中的 s ← χ_key。
        C++ 用 sampleHWT (汉明重量 h 的稀疏三元多项式)。
        本 demo 简化为系数 ∈ {-1, 0, 1} 的均匀采样。
        """
        return [self.rng.choice([-1, 0, 1]) for _ in range(self.par.N)]

    def sample_enc(self) -> Poly:
        """
        论文 Enc 中的 v ← χ_enc。
        C++ 用 sampleZO (zero-one 分布: 以概率 1/2 取 0，各 1/4 取 ±1)。
        本 demo 简化为 {-1, 0, 1} 均匀采样。
        """
        return [self.rng.choice([-1, 0, 1]) for _ in range(self.par.N)]

    def sample_error(self) -> Poly:
        """
        论文中的 e ← χ_err (误差分布)。
        C++ 用 sampleGauss (离散高斯，σ = 3.2)。
        本 demo 简化为 {-1, 0, 1} 均匀采样。

        注意: 极小的误差分布意味着密文噪声很小，
        在小参数下解密结果更容易正确，但这不代表安全性。
        """
        return [self.rng.choice([-1, 0, 1]) for _ in range(self.par.N)]

    def sample_uniform_rns(self, moduli: Sequence[int]) -> RNSPoly:
        """
        论文中 a ← R_Q 的均匀采样。
        在每个 R_{q_j} 上独立均匀采样多项式系数。
        """
        return {q: [self.rng.randrange(q) for _ in range(self.par.N)] for q in moduli}

    # ======================== 密钥生成: 论文 4.1 ========================

    def keygen(self) -> Tuple[SecretKey, PublicKey, EvalKey, EvalKey]:
        """
        论文第 4 节 KeyGen:
            1. s ← χ_key, sk = (1, s)
            2. pk = (b, a): a ← R_Q, e ← χ_err, b = -a*s + e
            3. evk_mult = KSGen(s^2, s)     // 乘法重线性化密钥
            4. evk_conj = KSGen(s_conj, s)  // 共轭密钥

        本 demo 同时生成旋转密钥。

        C++ 对应: Scheme::addEncKey, Scheme::addMultKey, Scheme::addConjKey
        """
        # Step 1: 生成 secret key
        s = self.sample_secret()
        sk = SecretKey(s=s)

        # Step 2: 生成 public key (论文 Enc 的公钥加密版本)
        # pk = (b, a), b = -a*s + e (mod Q_L)
        C_L = self.par.C(self.par.L)
        a_pk = self.sample_uniform_rns(C_L)
        e_pk = rns_from_poly(self.sample_error(), C_L)
        s_C = rns_from_poly(s, C_L)
        # b = -a*s + e
        minus_as = rns_scalar_mul(rns_mul(a_pk, s_C), -1)
        b_pk = rns_add(minus_as, e_pk)
        pk = PublicKey(b=b_pk, a=a_pk)

        # Step 3: 生成乘法评估密钥 evk_mult = KSGen(s^2, s)
        evk_mult = self.ksgen_mult(sk)

        # Step 4: 生成共轭密钥 evk_conj = KSGen(conjugate(s), s)
        evk_conj = self.ksgen_conjugate(sk)

        return sk, pk, evk_mult, evk_conj

    def ksgen_mult(self, sk: SecretKey) -> EvalKey:
        """
        论文第 4 节 KSGen(s^2, s): 乘法重线性化密钥。

        目标: 生成 evk = (b, a) ∈ R_{P*Q}^2，满足:
            b + a*s ≡ P*s^2 + e  (mod P*Q)

        这样在乘法 Mult 中，u2 * evk 解密后 ≈ P * u2 * s^2，
        再经 Algorithm 2 ModDown 近似除以 P，得到 u2 * s^2。

        C++ 对应: Scheme::addMultKey
            - sxsx = s*s (在 NTT 域)
            - evalAndEqual(sxsx): 对 Q 分量乘 P (论文中的 P 因子注入)
            - b = -a*s + (P*s^2 + e)

        本 demo 直接写 P*s^2，与 C++ 的 evalAndEqual 效果等价。
        """
        D_L = self.par.D(self.par.L)
        P = product(self.par.B)

        # a ← R_{D_L} (均匀随机)
        a = self.sample_uniform_rns(D_L)
        # s 和 s^2 在 D_L 基上的 RNS 表示
        s_D = rns_from_poly(sk.s, D_L)
        # s^2 先在整数意义下做环乘（大模数避免溢出），再转 RNS
        s2 = poly_mul(sk.s, sk.s, 10**18)
        s2_D = rns_from_poly(s2, D_L)
        # e ← χ_err
        e_D = rns_from_poly(self.sample_error(), D_L)

        # b = -a*s + P*s^2 + e  (mod each modulus in D_L)
        # 对比 C++: b = -a*s + evalAndEqual(s^2) + e，其中 evalAndEqual 等效于乘 P
        b = rns_add(
            rns_scalar_mul(rns_mul(a, s_D), -1),
            rns_add(rns_scalar_mul(s2_D, P), e_D)
        )
        return EvalKey(b=b, a=a)

    def ksgen_conjugate(self, sk: SecretKey) -> EvalKey:
        """
        论文第 4 节 KSGen(s_conj, s): 共轭密钥生成。

        目标: 生成 evk_conj = (b, a) ∈ R_{P*Q}^2，满足:
            b + a*s ≡ P * conjugate(s) + e  (mod P*Q)

        其中 conjugate(s) = s(X^{-1}) = s(X^{2N-1})。

        共轭操作将每个 slot 的值取复共轭。
        由于 key switching 需要知道 conjugate(s) 到 s 的转换关系，
        因此需要生成此评估密钥。

        C++ 对应: Scheme::addConjKey
            - sxconj = conjugate(s.sx)
            - evalAndEqual(sxconj): 乘 P
            - b = -a*s + (P*s_conj + e)
        """
        D_L = self.par.D(self.par.L)
        P = product(self.par.B)

        # a ← R_{D_L}
        a = self.sample_uniform_rns(D_L)
        s_D = rns_from_poly(sk.s, D_L)

        # conjugate(s) = s(X^{-1})
        s_conj = poly_conjugate(sk.s)
        s_conj_D = rns_from_poly(s_conj, D_L)

        # e ← χ_err
        e_D = rns_from_poly(self.sample_error(), D_L)

        # b = -a*s + P*conjugate(s) + e  (mod each modulus in D_L)
        b = rns_add(
            rns_scalar_mul(rns_mul(a, s_D), -1),
            rns_add(rns_scalar_mul(s_conj_D, P), e_D)
        )
        return EvalKey(b=b, a=a)

    def ksgen_rotate(self, sk: SecretKey, rot_slots: int) -> EvalKey:
        """
        论文第 4 节 KSGen(s_rot, s): 旋转密钥生成。

        目标: 生成 evk_rot = (b, a) ∈ R_{P*Q}^2，满足:
            b + a*s ≡ P * rotate(s, rotSlots) + e  (mod P*Q)

        其中 rotate(s) 是将 s(X) 映射为 s(X^{5^r}) 的操作
        (r = rotSlots, 5 是 Z_{2N}^* 的生成元)。

        C++ 对应: Scheme::addLeftRotKey
            - sxrot = leftRot(s.sx, rot)
            - evalAndEqual(sxrot): 乘 P
            - b = -a*s + (P*s_rot + e)
        """
        D_L = self.par.D(self.par.L)
        P = product(self.par.B)
        M = 2 * self.par.N

        # 计算旋转幂次: power = 5^rotSlots mod 2N
        # C++ Context 构造函数中 rotGroup[i] = 5^i mod M
        power = pow(5, rot_slots, M)

        # a ← R_{D_L}
        a = self.sample_uniform_rns(D_L)
        s_D = rns_from_poly(sk.s, D_L)

        # rotate(s): s(X) -> s(X^power)
        s_rot = poly_rotate(sk.s, power)
        s_rot_D = rns_from_poly(s_rot, D_L)

        # e ← χ_err
        e_D = rns_from_poly(self.sample_error(), D_L)

        # b = -a*s + P*rotate(s) + e
        b = rns_add(
            rns_scalar_mul(rns_mul(a, s_D), -1),
            rns_add(rns_scalar_mul(s_rot_D, P), e_D)
        )
        return EvalKey(b=b, a=a)

    # ======================== 编码/加密/解密: 论文 4.2 ========================

    def encrypt(self, values: Sequence[complex], pk: PublicKey,
                level: int | None = None) -> Ciphertext:
        """
        论文第 4 节 Enc_{pk}(m):

            v ← χ_enc, e0, e1 ← χ_err
            ct = v * pk + (m + e0, e1)
               = (v*b + m + e0,  v*a + e1)  ∈ R_{Q_l}^2

        编码 m 来自论文 2.1: m ≈ scale * σ^{-1}(values)。

        解密验证:
            c0 + c1*s = v*b + m + e0 + (v*a + e1)*s
                      = v*(b + a*s) + m + e0 + e1*s
                      ≈ v*e_pk + m + e0 + e1*s
                      ≈ m  (当噪声足够小时)

        C++ 对应: Scheme::encryptMsg(Plaintext)
        """
        if level is None:
            level = self.par.L
        C = self.par.C(level)

        # 论文 2.1: 编码 + 缩放
        m_poly = encode(values, self.par.N, self.par.scale)
        m = rns_from_poly(m_poly, C)

        # 论文 Enc: v, e0, e1 采样
        v = rns_from_poly(self.sample_enc(), C)
        e0 = rns_from_poly(self.sample_error(), C)
        e1 = rns_from_poly(self.sample_error(), C)

        # 限制 pk 到当前 level 的模数集
        pk_b = rns_restrict(pk.b, C)
        pk_a = rns_restrict(pk.a, C)

        # c0 = v*b + m + e0,  c1 = v*a + e1
        c0 = rns_add(rns_add(rns_mul(v, pk_b), m), e0)
        c1 = rns_add(rns_mul(v, pk_a), e1)
        return Ciphertext(c0=c0, c1=c1, level=level, slots=len(values))

    def encrypt_sk(self, values: Sequence[complex], sk: SecretKey,
                   level: int | None = None) -> Ciphertext:
        """
        论文第 4 节 秘密密钥加密 (可选):
            v ← χ_enc, e0 ← χ_err
            ct = (e0 + m, v)  使得 c0 + c1*s = e0 + m + v*s ≈ m

        注意: 本 demo 主要使用公钥加密，此函数仅供对比参考。
        """
        if level is None:
            level = self.par.L
        C = self.par.C(level)

        m_poly = encode(values, self.par.N, self.par.scale)
        m = rns_from_poly(m_poly, C)
        v = rns_from_poly(self.sample_enc(), C)
        e0 = rns_from_poly(self.sample_error(), C)

        c0 = rns_add(m, e0)
        c1 = v  # 简单地将 v 作为 c1
        return Ciphertext(c0=c0, c1=c1, level=level, slots=len(values))

    def decrypt_to_poly(self, ct: Ciphertext, sk: SecretKey) -> Poly:
        """
        论文第 4 节 Dec_{sk}(ct):

            对 ct = (c0, c1)，输出:
                m_noisy = c0 + c1 * s  (mod q_0)

        取最低层模数 q_0 上的中心代表元，得到带噪声的缩放明文。

        C++ 对应: Scheme::decryptMsg
            context.mul(mx, cipher.ax, secretKey.sx, 1);
            context.addAndEqual(mx, cipher.bx, 1);
        (注意 C++ 用 (ax, bx) = (c1, c0) 的顺序)
        """
        q0 = self.par.q[0]
        s_q0 = poly_mod(sk.s, q0)
        # m_noisy = c0 + c1*s (mod q0)
        raw = poly_add(ct.c0[q0], poly_mul(ct.c1[q0], s_q0, q0), q0)
        # 取中心代表元
        return [centered(x, q0) for x in raw]

    def decrypt(self, ct: Ciphertext, sk: SecretKey) -> List[complex]:
        """
        论文第 4 节: 先 Dec 再 Decode。

        Dec 得到带噪声的多项式 m_noisy，
        Decode (论文 2.1) 将其在 N/2 个 slot 上求值并除以 scale。

        C++ 对应: Scheme::decrypt -> decryptMsg + decode
        """
        return decode(self.decrypt_to_poly(ct, sk), self.par.N, self.par.scale)

    # ======================== 同态运算: 论文 4.3-4.4 ========================

    def add(self, lhs: Ciphertext, rhs: Ciphertext) -> Ciphertext:
        """
        论文第 4 节 Add(ct, ct'):

            ct_add = (c0 + d0, c1 + d1)

        解密: (c0+d0) + (c1+d1)*s = (c0+c1*s) + (d0+d1*s) ≈ m1 + m2

        每个 RNS 分量上逐项相加，不改变 level 和 scale。

        C++ 对应: Scheme::add
        """
        if lhs.level != rhs.level:
            raise ValueError("加法要求两个密文 level 相同")
        return Ciphertext(
            c0=rns_add(lhs.c0, rhs.c0),
            c1=rns_add(lhs.c1, rhs.c1),
            level=lhs.level,
            slots=lhs.slots,
        )

    def sub(self, lhs: Ciphertext, rhs: Ciphertext) -> Ciphertext:
        """
        论文第 4 节 Sub(ct, ct'):

            ct_sub = (c0 - d0, c1 - d1)

        解密: (c0-d0) + (c1-d1)*s = (c0+c1*s) - (d0+d1*s) ≈ m1 - m2

        C++ 对应: Scheme::sub
        """
        if lhs.level != rhs.level:
            raise ValueError("减法要求两个密文 level 相同")
        return Ciphertext(
            c0=rns_sub(lhs.c0, rhs.c0),
            c1=rns_sub(lhs.c1, rhs.c1),
            level=lhs.level,
            slots=lhs.slots,
        )

    def negate(self, ct: Ciphertext) -> Ciphertext:
        """
        论文第 4 节 Negate(ct):

            ct_neg = (-c0, -c1)

        解密: -c0 + (-c1)*s = -(c0+c1*s) ≈ -m

        C++ 对应: Scheme::negate
        """
        return Ciphertext(
            c0=rns_neg(ct.c0),
            c1=rns_neg(ct.c1),
            level=ct.level,
            slots=ct.slots,
        )

    def add_const(self, ct: Ciphertext, const: complex) -> Ciphertext:
        """
        论文第 4 节 AddConst(ct, c):

            将常数 c 编码为明文多项式，加到 c0 上。
            ct_new = (c0 + encode(c), c1)

        对于实数常数，只需在 c0 的常数项(0次项)加上 round(c * scale)。
        对于复数常数，还需在 X^{N/2} 项加上 round(imag(c) * scale)。

        C++ 对应: Scheme::addConst / Scheme::addConstAndEqual
            在 C++ 中，addConst 直接加 cnst*p 到 bx 的所有 NTT 系数上
            (NTT 域加常数 = 系数域加到常数项)
        """
        # 编码常数: 所有 slot 填入同一个值
        values = [const] * ct.slots
        m_const = encode(values, self.par.N, self.par.scale)
        m_rns = rns_from_poly(m_const, self.par.C(ct.level))

        return Ciphertext(
            c0=rns_add(ct.c0, m_rns),
            c1=ct.c1,  # c1 不变
            level=ct.level,
            slots=ct.slots,
        )

    def mult_by_const(self, ct: Ciphertext, const: complex) -> Ciphertext:
        """
        论文第 4 节 MultByConst(ct, c):

            将常数 c 编码为明文多项式，与密文逐分量相乘。
            ct_new = (c0 * encode(c), c1 * encode(c))

        解密: c0*m_c + c1*m_c*s = m_c*(c0+c1*s) ≈ c*m
        scale 从 Δ 变为 Δ^2 (需要 rescale)。

        C++ 对应: Scheme::multByConst
            tmpr = cnst * p, 然后 mulConst
        """
        values = [const] * ct.slots
        m_const = encode(values, self.par.N, self.par.scale)
        m_rns = rns_from_poly(m_const, self.par.C(ct.level))

        return Ciphertext(
            c0=rns_mul(ct.c0, m_rns),
            c1=rns_mul(ct.c1, m_rns),
            level=ct.level,
            slots=ct.slots,
        )

    def multiply(self, lhs: Ciphertext, rhs: Ciphertext,
                 evk: EvalKey) -> Ciphertext:
        """
        论文第 4 节 Mult_{evk}(ct, ct') 的完整步骤。

        设 ct = (c0, c1), ct' = (d0, d1)，解密分别是 c0+c1*s 与 d0+d1*s。
        直接相乘得到三项式:
            u0 = c0 * d0              (对应 secret 1)
            u1 = c0*d1 + c1*d0        (对应 secret s)
            u2 = c1 * d1              (对应 secret s^2)
        它对应扩展密钥 (1, s, s^2)。重线性化(relinearization)把 u2*s^2 转回二元密文。

        论文 Mult 的 5 个步骤:
            Step 1: 在 C 基上计算 u0, u1, u2
            Step 2: Algorithm 1 ModUp_{C->D}(u2)，扩展 u2 到 D = B ∪ C
            Step 3: u2_D * evk，在 D 基上与评估密钥相乘
            Step 4: Algorithm 2 ModDown_{D->C}，近似除以 P
            Step 5: 加回 u0, u1

        C++ 对应: Scheme::mult
            (C++ 用 Karatsuba: u1 = (c0+c1)(d0+d1) - u0 - u2，数学等价)

        注意: 乘法后 scale 从 Δ^2 (近似) 变为需要 rescale 回 Δ。
        """
        if lhs.level != rhs.level:
            raise ValueError("乘法要求两个密文 level 相同")

        level = lhs.level
        C = self.par.C(level)
        B = self.par.B

        # --- Step 1: 在 C 基上计算三项乘积 ---
        # u0 = c0 * d0
        u0 = rns_mul(lhs.c0, rhs.c0)
        # u1 = c0*d1 + c1*d0 (交叉项)
        u1 = rns_add(rns_mul(lhs.c0, rhs.c1), rns_mul(lhs.c1, rhs.c0))
        # u2 = c1 * d1
        u2 = rns_mul(lhs.c1, rhs.c1)

        # --- Step 2: Algorithm 1 ModUp_{C->D}(u2) ---
        # 给 u2 增加特殊基 B 分量，以便与 evk (存在于 D 基) 相乘
        u2_D = mod_up_C_to_D(u2, C, B)

        # --- Step 3: u2_D * evk (key switching 核心) ---
        # evk 满足 b + a*s ≈ P*s^2 (mod P*Q)
        # u2_D * evk.b ≈ u2 * P * s^2 (的 c0 贡献)
        # u2_D * evk.a ≈ u2 * P * s^2 (的 c1 贡献)
        D = B + C
        evk_b = rns_restrict(evk.b, D)
        evk_a = rns_restrict(evk.a, D)
        relin0_D = rns_mul(u2_D, evk_b)
        relin1_D = rns_mul(u2_D, evk_a)

        # --- Step 4: Algorithm 2 ModDown_{D->C} ---
        # 近似除以 P，将 P*s^2 的贡献还原为 s^2 的贡献
        relin0_C = mod_down_D_to_C(relin0_D, C, B)
        relin1_C = mod_down_D_to_C(relin1_D, C, B)

        # --- Step 5: 加回 u0, u1 ---
        return Ciphertext(
            c0=rns_add(u0, relin0_C),
            c1=rns_add(u1, relin1_C),
            level=level,
            slots=lhs.slots,
        )

    def rescale(self, ct: Ciphertext) -> Ciphertext:
        """
        论文第 4 节 RS_{l,l-1}(ct): Rescaling (重缩放/模切换)。

        输入 level=l 的密文，模数基 C_l = {q_0, ..., q_l}。
        丢弃最后一个模数 q_l，并对每个剩余 q_j 计算:

            c_i^{(j)} ← q_l^{-1} * (c_i^{(j)} - [c_i^{(l)}]_{q_j}) mod q_j

        这等价于:
            1. 从密文中减去 c_i mod q_l (使结果可被 q_l 整除)
            2. 除以 q_l
            3. 在剩余模数上取余

        效果:
            - 密文 scale 从 Δ^2 变回 ≈ Δ (乘法后必须 rescale)
            - level 减少 1
            - 引入近似误差: |c_i mod q_l| / q_l (论文 3.1 approximate basis 的误差源)

        C++ 对应: Context::reScaleAndEqual
            1. INTT(al) (将 NTT 表示转回系数)
            2. 对每个 q_i: rai = al mod q_i
            3. rai = (ai - rai) * q_l^{-1} mod q_i

        注意: 由于 q_l 只是"近似等于"固定 base q (论文 3.1)，
        rescale 后 scale 并非精确为 Δ，而是 Δ * (Δ / q_l)。
        当所有 q_j 相近时，这个比率接近 1，误差可接受。
        """
        if ct.level == 0:
            raise ValueError("level 0 的密文不能继续 rescale")

        q_drop = self.par.q[ct.level]  # 要丢弃的最高层模数 q_l
        C_new = self.par.C(ct.level - 1)  # 新的模数集 {q_0, ..., q_{l-1}}

        def rescale_component(comp: RNSPoly) -> RNSPoly:
            # 取出被丢弃模数上的分量 c^{(l)}
            dropped = comp[q_drop]
            out = {}
            for q in C_new:
                inv = inv_mod(q_drop, q)
                # c'^{(j)} = q_l^{-1} * (c^{(j)} - [c^{(l)}]_{q_j}) mod q_j
                out[q] = [((comp[q][i] - (dropped[i] % q)) * inv) % q
                          for i in range(self.par.N)]
            return out

        return Ciphertext(
            c0=rescale_component(ct.c0),
            c1=rescale_component(ct.c1),
            level=ct.level - 1,
            slots=ct.slots,
        )

    # ======================== 槽位操作: 共轭与旋转 ========================

    def conjugate(self, ct: Ciphertext, evk_conj: EvalKey) -> Ciphertext:
        """
        论文第 4 节 Conjugate(ct):

        目标: 将每个 slot 的值取复共轭。

        在多项式层面，共轭对应:
            a(X) -> a(X^{-1})

        对于密文 (c0, c1)，共轭后:
            ct_conj = (conjugate(c0), conjugate(c1))
        但 conjugate(c1) 的 secret 从 s 变成了 conjugate(s)，
        需要通过 key switching 将 conjugate(s) 转回 s。

        步骤:
            1. c0_conj = conjugate(c0),  c1_conj = conjugate(c1)
            2. 对 c1_conj 做 key switching (使用 evk_conj):
               - ModUp: c1_conj 从 C 扩展到 D
               - 乘 evk_conj: 得到 (c1_conj * evk.b, c1_conj * evk.a)
               - ModDown: 从 D 缩回 C
            3. 结果: (c0_conj + relin_0, relin_1)

        C++ 对应: Scheme::conjugate
        """
        level = ct.level
        C = self.par.C(level)
        B = self.par.B
        D = B + C

        # Step 1: 对 c0, c1 分别做共轭 (多项式 X -> X^{-1})
        c0_conj = {q: poly_conjugate(ct.c0[q]) for q in C}
        # 注意: 共轭后系数可能为负，需要取模
        c0_conj = {q: poly_mod(c0_conj[q], q) for q in C}

        c1_conj = {q: poly_conjugate(ct.c1[q]) for q in C}
        c1_conj = {q: poly_mod(c1_conj[q], q) for q in C}

        # Step 2: Key switching c1_conj (从 conjugate(s) 转回 s)
        # ModUp: 扩展到 D 基
        c1_conj_D = mod_up_C_to_D(c1_conj, C, B)

        # 乘评估密钥
        evk_b = rns_restrict(evk_conj.b, D)
        evk_a = rns_restrict(evk_conj.a, D)
        relin0_D = rns_mul(c1_conj_D, evk_b)
        relin1_D = rns_mul(c1_conj_D, evk_a)

        # ModDown: 缩回 C 基
        relin0_C = mod_down_D_to_C(relin0_D, C, B)
        relin1_C = mod_down_D_to_C(relin1_D, C, B)

        # Step 3: 组合结果
        return Ciphertext(
            c0=rns_add(c0_conj, relin0_C),
            c1=relin1_C,
            level=level,
            slots=ct.slots,
        )

    def left_rotate(self, ct: Ciphertext, rot_slots: int,
                    evk_rot: EvalKey) -> Ciphertext:
        """
        论文第 4 节 LeftRotate(ct, r):

        目标: 将 N/2 个 slot 循环左移 r 个位置。
        例如 [z0, z1, z2, z3] 左移 1 -> [z1, z2, z3, z0]。

        在多项式层面，旋转对应:
            a(X) -> a(X^{5^r})
        其中 5 是 Z_{2N}^* 的一个生成元。

        步骤 (类似 conjugate):
            1. c0_rot = rotate(c0),  c1_rot = rotate(c1)
            2. 对 c1_rot 做 key switching (使用 evk_rot)
            3. 结果: (c0_rot + relin_0, relin_1)

        C++ 对应: Scheme::leftRotateFast
        """
        level = ct.level
        C = self.par.C(level)
        B = self.par.B
        D = B + C
        M = 2 * self.par.N

        # 计算旋转幂次
        power = pow(5, rot_slots, M)

        # Step 1: 对 c0, c1 分别做旋转 (多项式 X -> X^power)
        c0_rot = {q: poly_mod(poly_rotate(ct.c0[q], power), q) for q in C}
        c1_rot = {q: poly_mod(poly_rotate(ct.c1[q], power), q) for q in C}

        # Step 2: Key switching c1_rot
        c1_rot_D = mod_up_C_to_D(c1_rot, C, B)

        evk_b = rns_restrict(evk_rot.b, D)
        evk_a = rns_restrict(evk_rot.a, D)
        relin0_D = rns_mul(c1_rot_D, evk_b)
        relin1_D = rns_mul(c1_rot_D, evk_a)

        relin0_C = mod_down_D_to_C(relin0_D, C, B)
        relin1_C = mod_down_D_to_C(relin1_D, C, B)

        # Step 3: 组合
        return Ciphertext(
            c0=rns_add(c0_rot, relin0_C),
            c1=relin1_C,
            level=level,
            slots=ct.slots,
        )


# ============================================================================
# 7. 辅助打印函数
# ============================================================================


def pretty_complex_vec(xs: Sequence[complex], digits: int = 4) -> List[str]:
    """将复数向量格式化为可读字符串列表。"""
    return [f"{z.real:.{digits}f}{z.imag:+.{digits}f}j" for z in xs]


def print_comparison(label: str, actual: List[complex],
                     expected: List[complex], digits: int = 4) -> None:
    """打印解密结果与期望值的对比。"""
    print(f"  {label}:")
    print(f"    解密结果: {pretty_complex_vec(actual, digits)}")
    print(f"    期望结果: {pretty_complex_vec(expected, digits)}")
    # 计算每个 slot 的误差
    errors = [abs(a - e) for a, e in zip(actual, expected)]
    max_err = max(errors)
    print(f"    最大误差: {max_err:.6f}")


# ============================================================================
# 8. 演示主流程
# ============================================================================


def main() -> None:
    """
    完整的 FullRNS-CKKS 教学演示。

    演示内容:
        1. Setup + KeyGen (含乘法密钥、共轭密钥、旋转密钥)
        2. 编码并加密两个 4-slot 复向量
        3. 同态加法 (Add)
        4. 同态减法 (Sub)
        5. 同态乘法 (Mult) + 重缩放 (Rescale)
        6. 常数加法 (AddConst)
        7. 常数乘法 (MultByConst) + Rescale
        8. 复共轭 (Conjugate)
        9. 槽位左旋转 (LeftRotate)

    注意:
        - 参数极小 (N=8, slots=4)，完全没有安全性
        - 乘法和 key switching 后误差可能很大
        - 这是预期行为，不是 bug
    """

    print("=" * 70)
    print("FullRNS-HEAAN / CKKS 教学演示")
    print("论文: A Full RNS Variant of Approximate Homomorphic Encryption")
    print("=" * 70)

    # ---- Step 1: Setup + KeyGen ----
    print("\n[Step 1] Setup + KeyGen (论文第 4 节)")
    par = Params()
    he = FullRNSCKKSDemo(par)
    print(f"  N = {par.N}, slots = {par.N // 2}, L = {par.L}")
    print(f"  scale (Δ) = {par.scale} = 2^{int(math.log2(par.scale))}")
    print(f"  密文模数 q = {par.q}")
    print(f"  特殊模数 B = {par.B}")
    print(f"  P = prod(B) = {product(par.B)}")

    sk, pk, evk_mult, evk_conj = he.keygen()
    print(f"  secret key s = {sk.s}")

    # 生成旋转密钥 (rot_slots = 1)
    evk_rot1 = he.ksgen_rotate(sk, rot_slots=1)
    print("  已生成: pk, evk_mult, evk_conj, evk_rot(1)")

    # ---- 测试数据 ----
    z = [1.25 + 0.50j, -0.75 + 0.25j, 0.50 - 0.50j, -1.00 - 0.25j]
    w = [0.50 - 0.25j,  1.00 + 0.00j, -0.25 + 0.75j,  0.75 + 0.50j]

    print(f"\n  输入向量 z = {pretty_complex_vec(z)}")
    print(f"  输入向量 w = {pretty_complex_vec(w)}")

    # ---- Step 2: 加密 ----
    print("\n[Step 2] Enc (论文第 4 节)")
    ct_z = he.encrypt(z, pk)
    ct_w = he.encrypt(w, pk)
    print(f"  ct_z: level={ct_z.level}, slots={ct_z.slots}")
    print(f"  ct_w: level={ct_w.level}, slots={ct_w.slots}")

    # 验证单独解密
    dec_z = he.decrypt(ct_z, sk)
    dec_w = he.decrypt(ct_w, sk)
    print_comparison("Dec(ct_z)", dec_z, z)
    print_comparison("Dec(ct_w)", dec_w, w)

    # ---- Step 3: 同态加法 ----
    print("\n[Step 3] Add (论文第 4 节)")
    ct_add = he.add(ct_z, ct_w)
    dec_add = he.decrypt(ct_add, sk)
    expected_add = [a + b for a, b in zip(z, w)]
    print_comparison("Dec(ct_z + ct_w)", dec_add, expected_add)

    # ---- Step 4: 同态减法 ----
    print("\n[Step 4] Sub (论文第 4 节)")
    ct_sub = he.sub(ct_z, ct_w)
    dec_sub = he.decrypt(ct_sub, sk)
    expected_sub = [a - b for a, b in zip(z, w)]
    print_comparison("Dec(ct_z - ct_w)", dec_sub, expected_sub)

    # ---- Step 5: 同态取负 ----
    print("\n[Step 5] Negate (论文第 4 节)")
    ct_neg = he.negate(ct_z)
    dec_neg = he.decrypt(ct_neg, sk)
    expected_neg = [-a for a in z]
    print_comparison("Dec(-ct_z)", dec_neg, expected_neg)

    # ---- Step 6: 常数加法 ----
    print("\n[Step 6] AddConst (论文第 4 节)")
    const_add = 0.5 + 0.25j
    ct_addc = he.add_const(ct_z, const_add)
    dec_addc = he.decrypt(ct_addc, sk)
    expected_addc = [a + const_add for a in z]
    print(f"  常数 = {const_add}")
    print_comparison("Dec(ct_z + 0.5+0.25j)", dec_addc, expected_addc)

    # ---- Step 7: 同态乘法 + Rescale ----
    print("\n[Step 7] Mult + Rescale (论文第 4 节)")
    ct_mul = he.multiply(ct_z, ct_w, evk_mult)
    print(f"  乘法后: level={ct_mul.level} (scale ≈ Δ^2 = {par.scale**2})")

    # Rescale: 丢弃最高层模数，scale 从 Δ^2 回到 ≈ Δ
    ct_mul_rs = he.rescale(ct_mul)
    print(f"  Rescale 后: level={ct_mul_rs.level} (scale ≈ Δ = {par.scale})")

    dec_mul = he.decrypt(ct_mul_rs, sk)
    expected_mul = [a * b for a, b in zip(z, w)]
    print_comparison("Dec(ct_z * ct_w)", dec_mul, expected_mul)
    print("  注意: 乘法涉及 ModUp + key switching + ModDown + Rescale，")
    print("  在小参数下误差会显著增大，这是正常的。")

    # ---- Step 8: 常数乘法 + Rescale ----
    print("\n[Step 8] MultByConst + Rescale (论文第 4 节)")
    const_mul = 2.0 + 0.5j
    ct_mulc = he.mult_by_const(ct_z, const_mul)
    ct_mulc_rs = he.rescale(ct_mulc)
    dec_mulc = he.decrypt(ct_mulc_rs, sk)
    expected_mulc = [a * const_mul for a in z]
    print(f"  常数 = {const_mul}")
    print_comparison("Dec(ct_z * (2+0.5j))", dec_mulc, expected_mulc)

    # ---- Step 9: 复共轭 ----
    print("\n[Step 9] Conjugate (论文第 4 节)")
    ct_conj = he.conjugate(ct_z, evk_conj)
    dec_conj = he.decrypt(ct_conj, sk)
    expected_conj = [a.conjugate() for a in z]
    print_comparison("Dec(conjugate(ct_z))", dec_conj, expected_conj)
    print("  共轭操作将每个 slot 取复共轭: a+bi -> a-bi")

    # ---- Step 10: 槽位左旋转 ----
    print("\n[Step 10] LeftRotate (论文第 4 节)")
    ct_rot = he.left_rotate(ct_z, rot_slots=1, evk_rot=evk_rot1)
    dec_rot = he.decrypt(ct_rot, sk)
    # 期望: 左移 1 位, [z0,z1,z2,z3] -> [z1,z2,z3,z0]
    expected_rot = z[1:] + z[:1]
    print_comparison("Dec(rotate(ct_z, 1))", dec_rot, expected_rot)
    print("  旋转操作将 slot 循环左移: [z0,z1,z2,z3] -> [z1,z2,z3,z0]")

    # ---- 总结 ----
    print("\n" + "=" * 70)
    print("演示完成！")
    print("=" * 70)
    print("\n关键观察:")
    print("  1. 加法/减法的误差很小(仅来自编码舍入和加密噪声)")
    print("  2. 乘法 + Rescale 的误差显著增大，因为:")
    print("     - ModUp/ModDown 引入 fast basis conversion 的近似误差")
    print("     - Rescale 丢弃 q_l 时引入 round-off 误差")
    print("     - 论文 3.1 的 approximate basis (各 q_j 不完全相等) 放大误差")
    print("  3. 共轭和旋转也涉及 key switching (ModUp + ModDown)，误差类似")
    print("  4. 在真实参数下 (N ≥ 2^13, q_j ≈ 2^60)，这些误差相对 Δ 可忽略")
    print("\n论文各节与本 demo 代码的对应关系详见文件头部注释和各函数文档字符串。")


if __name__ == "__main__":
    main()
