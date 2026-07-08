"""
FullRNS-HEAAN / CKKS 教学 demo
================================

这个文件不是为了安全性或性能，而是为了把论文
“A Full RNS Variant of Approximate Homomorphic Encryption”(2018/931)
中的核心方案步骤，用尽量小、尽量可读的 Python 代码串起来。

本仓库 C++ 代码的结构大致对应：
    - src/Context.*      : 环运算、NTT、RNS 基转换、raise/back/rescale
    - src/Scheme.*       : KeyGen、加解密、同态加法/乘法、重缩放
    - src/SecretKey.*    : 小汉明重量 secret key

本 demo 和论文/代码的主要对应关系：
    - 论文 2.1: CKKS/HEAAN 的近似数编码、缩放、同态加乘和 rescale 思想
    - 论文 2.2: RNS 表示 [a]_B
    - 论文 2.3: Fast Basis Conversion Conv_C->B / Conv_B->C
    - 论文 3.1: approximate basis，即 q_j 近似同一个固定 base q
    - 论文 3.2: Algorithm 1 ModUp、Algorithm 2 ModDown
    - 论文 4  : Setup、KSGen、KeyGen、Enc、Dec、Add、Mult、RS

注意：
    1. 真实库会把每个 RNS 分量再转为 NTT 表示以做快速多项式乘法；
       这里为了教学，把多项式保存在“系数表示”中，并逐模数做环乘法。
       这不改变 RNS 公式，只是省掉 NTT 优化层。
    2. 参数极小：N=8，因此可打包 slots=N/2=4 个复数。
       这完全没有安全性，误差也会明显偏大。
    3. 为了便于重复阅读，随机性使用固定 seed。
"""

from __future__ import annotations

from dataclasses import dataclass
import cmath
import random
from typing import Dict, Iterable, List, Sequence, Tuple


Poly = List[int]
RNSPoly = Dict[int, Poly]          # modulus -> polynomial coefficients modulo modulus
RNSTuple = Tuple[RNSPoly, RNSPoly] # ciphertext component pair (c0, c1)


# ---------------------------------------------------------------------------
# 1. 小工具：模逆、中心代表、多项式环 R = Z[X]/(X^N + 1)
# ---------------------------------------------------------------------------


def inv_mod(a: int, q: int) -> int:
    """返回 a^{-1} mod q。Python 3.8+ 支持 pow(a, -1, q)。"""
    return pow(a % q, -1, q)


def centered(x: int, q: int) -> int:
    """论文第 2 节记号 [a]_q: 取 (-q/2, q/2] 中的中心代表。"""
    x %= q
    return x if x <= q // 2 else x - q


def poly_zero(N: int) -> Poly:
    return [0] * N


def poly_add(a: Poly, b: Poly, q: int) -> Poly:
    return [(x + y) % q for x, y in zip(a, b)]


def poly_sub(a: Poly, b: Poly, q: int) -> Poly:
    return [(x - y) % q for x, y in zip(a, b)]


def poly_scalar_mul(a: Poly, s: int, q: int) -> Poly:
    return [(x * s) % q for x in a]


def poly_mul(a: Poly, b: Poly, q: int) -> Poly:
    """
    环 R_q = Z_q[X]/(X^N + 1) 中的乘法。

    对普通卷积中次数 >= N 的项，使用 X^N = -1 折回：
        X^{N+t} = -X^t.

    论文第 4 节所有密文分量都属于 R_{q_j} 或 R_{p_i}。
    """
    N = len(a)
    out = [0] * N
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            deg = i + j
            if deg < N:
                out[deg] += ai * bj
            else:
                out[deg - N] -= ai * bj
    return [x % q for x in out]


def poly_mod(poly: Poly, q: int) -> Poly:
    return [x % q for x in poly]


def rns_from_poly(poly: Poly, moduli: Sequence[int]) -> RNSPoly:
    """论文 2.2 的 RNS 表示 [a]_B，按系数对每个模数取余。"""
    return {q: poly_mod(poly, q) for q in moduli}


def rns_add(a: RNSPoly, b: RNSPoly) -> RNSPoly:
    return {q: poly_add(a[q], b[q], q) for q in a}


def rns_sub(a: RNSPoly, b: RNSPoly) -> RNSPoly:
    return {q: poly_sub(a[q], b[q], q) for q in a}


def rns_mul(a: RNSPoly, b: RNSPoly) -> RNSPoly:
    return {q: poly_mul(a[q], b[q], q) for q in a}


def rns_scalar_mul(a: RNSPoly, s: int) -> RNSPoly:
    return {q: poly_scalar_mul(a[q], s, q) for q in a}


def rns_restrict(a: RNSPoly, moduli: Sequence[int]) -> RNSPoly:
    return {q: a[q][:] for q in moduli}


# ---------------------------------------------------------------------------
# 2. 论文 2.3: Fast Basis Conversion
# ---------------------------------------------------------------------------


def product(xs: Iterable[int]) -> int:
    r = 1
    for x in xs:
        r *= x
    return r


def fast_convert(poly_in_basis: RNSPoly, src: Sequence[int], dst: Sequence[int]) -> RNSPoly:
    """
    论文 2.3 Fast Basis Conversion:

        Conv_{src -> dst}([a]_src)_i
          = sum_j ([a^{(j)} * qhat_j^{-1}]_{q_j} * qhat_j) mod p_i

    其中 qhat_j = prod_{t != j} q_t。

    重要性质：
        输出通常不是精确整数 a 在新基上的 residue，而是
        a + Q * e 的 residue；e 较小。
        论文 3.2 的 ModUp/ModDown 正是利用这个“近似但足够好”的性质。
    """
    N = len(next(iter(poly_in_basis.values())))
    src_prod = product(src)
    out = {p: [0] * N for p in dst}

    for j, qj in enumerate(src):
        qhat = src_prod // qj
        qhat_inv_mod_qj = inv_mod(qhat, qj)
        scaled_coeffs = [(c * qhat_inv_mod_qj) % qj for c in poly_in_basis[qj]]

        for p in dst:
            qhat_mod_p = qhat % p
            for n in range(N):
                out[p][n] = (out[p][n] + scaled_coeffs[n] * qhat_mod_p) % p
    return out


def mod_up_C_to_D(a_C: RNSPoly, C: Sequence[int], B: Sequence[int]) -> RNSPoly:
    """
    论文 3.2 Algorithm 1: Approximate Modulus Raising ModUp_{C->D}

    输入:  [a]_C, C={q_0,...,q_l}
    输出:  [a']_D, D=B union C, 且 a' = a (mod Q_l), |a'| 约小于 P*Q_l

    步骤:
        1. 用 Conv_{C->B} 生成特殊基 B 上的 residue。
        2. 拼接 B 分量和原 C 分量。
    """
    a_B = fast_convert(a_C, C, B)
    return {**a_B, **{q: a_C[q][:] for q in C}}


def mod_down_D_to_C(a_D: RNSPoly, C: Sequence[int], B: Sequence[int]) -> RNSPoly:
    """
    论文 3.2 Algorithm 2: Approximate Modulus Reduction ModDown_{D->C}

    目标:
        输入 [b]_D, D=B union C，输出 [b']_C，使 b' 约等于 P^{-1} * b。

    公式对应论文 Algorithm 2 第 2-5 行：
        a_C = Conv_{B->C}([b]_B)
        b'_j = P^{-1} * (b_{C,j} - a_{C,j}) mod q_j

    直觉:
        先从 B 分量恢复一个“小的 a”，满足 a = b (mod P)。
        于是 b-a 可被 P 整除；在 C 分量上乘 P^{-1} 就等价于除以 P。
    """
    P = product(B)
    a_B = {p: a_D[p] for p in B}
    approx_a_C = fast_convert(a_B, B, C)
    out = {}
    for q in C:
        pinv = inv_mod(P, q)
        out[q] = [((a_D[q][i] - approx_a_C[q][i]) * pinv) % q for i in range(len(a_D[q]))]
    return out


# ---------------------------------------------------------------------------
# 3. 论文 2.1: CKKS 编码/解码，打包 slots=N/2 个复数
# ---------------------------------------------------------------------------


def primitive_root_exponents_for_slots(N: int) -> List[int]:
    """
    论文 2.1 使用 canonical embedding 的一半坐标。

    论文文字写作 a(zeta), a(zeta^5), ...；
    对 N=8，取指数 [1,5,9,13]，正好是 2N 次单位根中一组共轭代表。
    """
    return [(1 + 4 * j) % (2 * N) for j in range(N // 2)]


def encode(values: Sequence[complex], N: int, scale: int) -> Poly:
    """
    论文 2.1 的 Encode = sigma^{-1} 后乘 scale 并舍入到 R。

    这里直接做插值：
        给定 slots=N/2 个复数 z_j，
        要找 m(X) 使 m(zeta^{e_j}) ~= scale * z_j。

    为保证系数接近实数，对未显式存储的共轭坐标填入 conjugate(z_j)。
    最后用逆 DFT 风格的公式恢复 degree < N 的系数。
    """
    slots = N // 2
    if len(values) != slots:
        raise ValueError(f"N={N} 时需要 {slots} 个 slot")

    zeta = cmath.exp(1j * cmath.pi / N)
    selected = primitive_root_exponents_for_slots(N)

    y_by_exp = {}
    for e, z in zip(selected, values):
        y_by_exp[e] = scale * z
        y_by_exp[(-e) % (2 * N)] = scale * z.conjugate()

    odd_exponents = [2 * r + 1 for r in range(N)]
    coeffs = []
    for k in range(N):
        ck = sum(y_by_exp[e] * (zeta ** (-e * k)) for e in odd_exponents) / N
        coeffs.append(int(round(ck.real)))
    return coeffs


def decode(poly: Poly, N: int, scale: int) -> List[complex]:
    """论文 2.1 的 Decode = sigma(m) / scale，只取 N/2 个 slot。"""
    zeta = cmath.exp(1j * cmath.pi / N)
    out = []
    for e in primitive_root_exponents_for_slots(N):
        val = sum(poly[k] * (zeta ** (e * k)) for k in range(N)) / scale
        out.append(val)
    return out


# ---------------------------------------------------------------------------
# 4. 参数、密钥、密文：论文第 4 节 Setup / KeyGen / Enc / Dec
# ---------------------------------------------------------------------------


@dataclass
class Params:
    """
    论文第 4 节 Setup(q,L,epsilon;lambda) 的教学参数。

    真实论文要求:
        q_j 约等于固定 base q，且 q_j = 1 mod 2N 以支持 NTT。

    本 demo:
        N=8，所以 2N=16；
        q_i 和 p_i 都选择 1 mod 16 的小素数。
        q_i 大小都在 2^15 附近，因此是 approximate basis 的玩具版。
    """

    N: int = 8
    L: int = 2
    scale: int = 2**8
    q: Tuple[int, ...] = (40961, 65537, 114689) # ciphertext moduli q_0,q_1,q_2
    B: Tuple[int, ...] = (147457, 163841)       # special primes p_0,p_1 for key switching

    def C(self, level: int) -> Tuple[int, ...]:
        return self.q[: level + 1]

    def D(self, level: int) -> Tuple[int, ...]:
        return self.B + self.C(level)


@dataclass
class SecretKey:
    s: Poly


@dataclass
class PublicKey:
    b: RNSPoly
    a: RNSPoly


@dataclass
class EvalKey:
    b: RNSPoly
    a: RNSPoly


@dataclass
class Ciphertext:
    c0: RNSPoly
    c1: RNSPoly
    level: int
    slots: int


class FullRNSCKKSDemo:
    def __init__(self, params: Params, seed: int = 7):
        self.par = params
        self.rng = random.Random(seed)

    # -------------------------- sampling distributions --------------------
    # 论文第 4 节 Setup 选择 key、enc、err 分布；附录 A 讨论其噪声。
    # 这里使用极小噪声，方便肉眼观察，不用于安全。

    def sample_secret(self) -> Poly:
        """玩具 secret key: 系数来自 {-1,0,1}，类似小汉明重量 key。"""
        return [self.rng.choice([-1, 0, 1]) for _ in range(self.par.N)]

    def sample_enc(self) -> Poly:
        """论文 Enc 中的 v <- enc；这里也取 {-1,0,1}。"""
        return [self.rng.choice([-1, 0, 1]) for _ in range(self.par.N)]

    def sample_error(self) -> Poly:
        """论文 err 分布；这里取 {-1,0,1} 的微小误差。"""
        return [self.rng.choice([-1, 0, 1]) for _ in range(self.par.N)]

    def sample_uniform_rns(self, moduli: Sequence[int]) -> RNSPoly:
        """在每个 R_q 上采样均匀多项式。"""
        return {q: [self.rng.randrange(q) for _ in range(self.par.N)] for q in moduli}

    # -------------------------- key generation ----------------------------

    def keygen(self) -> Tuple[SecretKey, PublicKey, EvalKey]:
        """
        论文第 4 节 KeyGen:

            s <- key, sk=(1,s)
            pk=(b,a), b = -a*s + e       (本代码用 c0 + c1*s 解密)
            evk = KSGen(s^2, s)

        论文符号中有时写成 <ct,sk>，本 demo 统一用：
            Dec(ct) = c0 + c1*s.
        """
        s = self.sample_secret()
        sk = SecretKey(s=s)

        C_L = self.par.C(self.par.L)
        a_pk = self.sample_uniform_rns(C_L)
        e_pk = rns_from_poly(self.sample_error(), C_L)
        s_C = rns_from_poly(s, C_L)
        minus_a_s = rns_scalar_mul(rns_mul(a_pk, s_C), -1)
        b_pk = rns_add(minus_a_s, e_pk)
        pk = PublicKey(b=b_pk, a=a_pk)

        evk = self.ksgen_s_square_to_s(sk)
        return sk, pk, evk

    def ksgen_s_square_to_s(self, sk: SecretKey) -> EvalKey:
        """
        论文第 4 节 KSGen(s1,s2)，乘法重线性化使用 s1=s^2, s2=s。

        目标是发布 evk=(b,a) in R_{P*Q}^2，使：
            b + a*s = P*s^2 + e    (mod P*Q)

        这样在乘法中 d2*evk 解密到 P*d2*s^2；
        随后 Algorithm 2 ModDown 近似除以 P，得到 d2*s^2 的普通密文贡献。
        """
        D_L = self.par.D(self.par.L)
        P = product(self.par.B)

        a = self.sample_uniform_rns(D_L)
        s_D = rns_from_poly(sk.s, D_L)
        s2 = poly_mul(sk.s, sk.s, 10**18) # 先在整数意义下得到环乘结构
        s2_D = rns_from_poly(s2, D_L)
        e_D = rns_from_poly(self.sample_error(), D_L)

        # b = -a*s + P*s^2 + e  (mod each modulus in D)
        b = rns_add(rns_scalar_mul(rns_mul(a, s_D), -1), rns_scalar_mul(s2_D, P))
        b = rns_add(b, e_D)
        return EvalKey(b=b, a=a)

    # -------------------------- encode/encrypt/decrypt --------------------

    def encrypt(self, values: Sequence[complex], pk: PublicKey, level: int | None = None) -> Ciphertext:
        """
        论文第 4 节 Encpk(m):

            v <- enc, e0,e1 <- err
            ct = v*pk + (m+e0, e1)  in prod_j R_{q_j}^2

        编码 m 来自论文 2.1：m ~= scale * sigma^{-1}(values)。
        """
        if level is None:
            level = self.par.L
        C = self.par.C(level)

        m_poly = encode(values, self.par.N, self.par.scale)
        m = rns_from_poly(m_poly, C)
        v = rns_from_poly(self.sample_enc(), C)
        e0 = rns_from_poly(self.sample_error(), C)
        e1 = rns_from_poly(self.sample_error(), C)

        pk_b = rns_restrict(pk.b, C)
        pk_a = rns_restrict(pk.a, C)
        c0 = rns_add(rns_add(rns_mul(v, pk_b), m), e0)
        c1 = rns_add(rns_mul(v, pk_a), e1)
        return Ciphertext(c0=c0, c1=c1, level=level, slots=len(values))

    def decrypt_to_poly(self, ct: Ciphertext, sk: SecretKey) -> Poly:
        """
        论文第 4 节 Decsk(ct):

            对 ct=(c0,c1)，输出 c0 + c1*s (mod q0)。

        CKKS 解码只需要低层 q0 上的中心代表；真实库中会先 INTT。
        """
        q0 = self.par.q[0]
        s_q0 = poly_mod(sk.s, q0)
        raw = poly_add(ct.c0[q0], poly_mul(ct.c1[q0], s_q0, q0), q0)
        return [centered(x, q0) for x in raw]

    def decrypt(self, ct: Ciphertext, sk: SecretKey) -> List[complex]:
        return decode(self.decrypt_to_poly(ct, sk), self.par.N, self.par.scale)

    # -------------------------- homomorphic operations --------------------

    def add(self, lhs: Ciphertext, rhs: Ciphertext) -> Ciphertext:
        """
        论文第 4 节 Add:
            每个 RNS 分量上逐项相加。
        """
        if lhs.level != rhs.level:
            raise ValueError("这个教学 demo 只相加相同 level 的密文")
        return Ciphertext(
            c0=rns_add(lhs.c0, rhs.c0),
            c1=rns_add(lhs.c1, rhs.c1),
            level=lhs.level,
            slots=lhs.slots,
        )

    def multiply(self, lhs: Ciphertext, rhs: Ciphertext, evk: EvalKey) -> Ciphertext:
        """
        论文第 4 节 Multevk(ct,ct') 的 5 个步骤。

        设 ct=(c0,c1), ct'=(d0,d1)，解密分别是 c0+c1*s 与 d0+d1*s。
        直接相乘会得到三项：
            u0 = c0*d0
            u1 = c0*d1 + c1*d0
            u2 = c1*d1
        它对应 secret key (1,s,s^2)。重线性化把 u2*s^2 转回两项密文。
        """
        if lhs.level != rhs.level:
            raise ValueError("这个教学 demo 只相乘相同 level 的密文")

        level = lhs.level
        C = self.par.C(level)
        B = self.par.B

        # Step 1: 在当前 C 基上计算三项乘积。
        u0 = rns_mul(lhs.c0, rhs.c0)
        u1 = rns_add(rns_mul(lhs.c0, rhs.c1), rns_mul(lhs.c1, rhs.c0))
        u2 = rns_mul(lhs.c1, rhs.c1)

        # Step 2: Algorithm 1 ModUp_{C->D}(u2)，给 u2 增加特殊基 B 分量。
        u2_D = mod_up_C_to_D(u2, C, B)

        # Step 3: u2_D * evk，在 D=B union C 的每个模数上相乘。
        D = B + C
        evk_b = rns_restrict(evk.b, D)
        evk_a = rns_restrict(evk.a, D)
        relin0_D = rns_mul(u2_D, evk_b)
        relin1_D = rns_mul(u2_D, evk_a)

        # Step 4: Algorithm 2 ModDown_{D->C}，近似除以 P。
        relin0_C = mod_down_D_to_C(relin0_D, C, B)
        relin1_C = mod_down_D_to_C(relin1_D, C, B)

        # Step 5: 加回 u0,u1，得到普通二元密文。
        return Ciphertext(
            c0=rns_add(u0, relin0_C),
            c1=rns_add(u1, relin1_C),
            level=level,
            slots=lhs.slots,
        )

    def rescale(self, ct: Ciphertext) -> Ciphertext:
        """
        论文第 4 节 RS(ct):

        输入 level=l 的密文，模数基 C_l={q0,...,ql}。
        丢弃最后一个模数 q_l，并对每个剩余 q_j 计算：

            c_i^{(j)} <- q_l^{-1} * (c_i^{(j)} - c_i^{(l)}) mod q_j

        这等价于把加密明文近似除以 q_l。
        在 CKKS 固定点语义中，乘法后 scale 从 Delta^2 变成 Delta；
        本论文第 3.1 节允许 q_l 只是“近似等于”固定 base q，因此会额外引入
        approximate basis 的小误差。
        """
        if ct.level == 0:
            raise ValueError("level 0 不能继续 rescale")

        q_drop = self.par.q[ct.level]
        C_old = self.par.C(ct.level)
        C_new = self.par.C(ct.level - 1)

        def rescale_component(comp: RNSPoly) -> RNSPoly:
            dropped = comp[q_drop]
            out = {}
            for q in C_new:
                inv = inv_mod(q_drop, q)
                out[q] = [((comp[q][i] - (dropped[i] % q)) * inv) % q for i in range(self.par.N)]
            return out

        # C_old 变量保留在这里主要是为了和论文 C_l 记号对应，便于断点观察。
        _ = C_old
        return Ciphertext(
            c0=rescale_component(ct.c0),
            c1=rescale_component(ct.c1),
            level=ct.level - 1,
            slots=ct.slots,
        )


# ---------------------------------------------------------------------------
# 5. 一个极小演示流程
# ---------------------------------------------------------------------------


def pretty_complex_vec(xs: Sequence[complex], digits: int = 4) -> List[str]:
    return [f"{z.real:.{digits}f}{z.imag:+.{digits}f}j" for z in xs]


def main() -> None:
    """
    这个 main 是阅读入口：

        1. Setup/KeyGen
        2. 编码并加密两个 4-slot 复向量
        3. 同态加法
        4. 同态乘法
        5. rescale 一层
        6. 解密观察近似结果

    用户要求中说明“不需要运行和调试”，所以这里仅提供可执行入口。
    小参数下误差可能很大，尤其是乘法+重线性化+rescale 后。
    """
    par = Params()
    he = FullRNSCKKSDemo(par)

    sk, pk, evk = he.keygen()

    z = [1.25 + 0.50j, -0.75 + 0.25j, 0.50 - 0.50j, -1.00 - 0.25j]
    w = [0.50 - 0.25j, 1.00 + 0.00j, -0.25 + 0.75j, 0.75 + 0.50j]

    ct_z = he.encrypt(z, pk)
    ct_w = he.encrypt(w, pk)

    ct_add = he.add(ct_z, ct_w)
    dec_add = he.decrypt(ct_add, sk)

    ct_mul = he.multiply(ct_z, ct_w, evk)
    ct_mul_rs = he.rescale(ct_mul)
    dec_mul = he.decrypt(ct_mul_rs, sk)

    print("原始 z:       ", pretty_complex_vec(z))
    print("原始 w:       ", pretty_complex_vec(w))
    print("解密 z+w:     ", pretty_complex_vec(dec_add))
    print("期望 z+w:     ", pretty_complex_vec([a + b for a, b in zip(z, w)]))
    print("解密 z*w:     ", pretty_complex_vec(dec_mul))
    print("期望 z*w:     ", pretty_complex_vec([a * b for a, b in zip(z, w)]))


if __name__ == "__main__":
    main()
