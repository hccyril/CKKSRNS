# FullRNS-HEAAN 方案说明文档

> **论文**: "A Full RNS Variant of Approximate Homomorphic Encryption" (SAC 2018)  
> **论文链接**: https://eprint.iacr.org/2018/931  
> **作者**: Yongsoo Song, Yongsoo Hwang, Miran Kim, Yongdae Kim 等 (CryptoLab / SnuCL)  
> **仓库**: HEAAN 的 RNS 变体实现  
> **许可证**: CC BY-NC 3.0 (非商业用途)

---

## 目录

1. [项目总体架构](#1-项目总体架构)
2. [论文核心方案概述](#2-论文核心方案概述)
3. [项目结构详解](#3-项目结构详解)
4. [关键代码文件详解](#4-关键代码文件详解)
5. [代码与论文精准映射关系](#5-代码与论文精准映射关系)
6. [核心算法流程图解](#6-核心算法流程图解)
7. [数学基础速查](#7-数学基础速查)
8. [关键数据结构](#8-关键数据结构)

---

## 1. 项目总体架构

### 1.1 项目定位

本仓库是论文 **"A Full RNS Variant of Approximate Homomorphic Encryption"** 的 C++ 参考实现。该论文是原始 HEAAN (CKKS) 方案的 RNS 变体，核心贡献是：

- **消除大整数运算**: 原始 HEAAN 需要任意精度大整数运算；本方案使用 RNS (Residue Number System) 表示，所有运算在固定位宽的模数上完成
- **快速基转换**: 提出高效的 `Conv_{B→C}` 算法，替代耗时的 CRT 重建
- **近似基 (Approximate Basis)**: 允许所有密文模数 `q_j` 近似相等，简化参数选择
- **近似 ModUp/ModDown**: 利用快速基转换的"近似但足够好"性质，实现高效的模数提升与降低

### 1.2 编译与运行

```bash
# 编译 HE 库
cd lib && make clean && make all

# 编译并运行测试
cd run && make clean && make && ./FRNSHEAAN
```

**注意**: 本代码仅在 GCC 下编译通过，因为使用了 `unsigned __int128` 类型（GCC 特有扩展）。

### 1.3 架构分层

```
┌─────────────────────────────────────────────────┐
│                  应用层 (main.cpp)                │
├─────────────────────────────────────────────────┤
│          算法层 (SchemeAlgo)                      │
│    power / inverse / exponent / sigmoid / FFT    │
├─────────────────────────────────────────────────┤
│          方案层 (Scheme)                          │
│  Enc / Dec / Add / Mult / RS / Rotate / Conj    │
├─────────────────────────────────────────────────┤
│          环运算层 (Context)                       │
│  NTT / INTT / encode / decode / 采样 / RNS ops  │
├─────────────────────────────────────────────────┤
│          基础层 (Numb / EvaluatorUtils)          │
│  模运算 / Barrett / 素数检测 / 随机数 / 旋转     │
├─────────────────────────────────────────────────┤
│          数据层 (Ciphertext / Plaintext / Key)    │
│  密文结构 / 明文结构 / 密钥结构                   │
└─────────────────────────────────────────────────┘
```

---

## 2. 论文核心方案概述

### 2.1 论文解决的问题

原始 CKKS/HEAAN 方案的瓶颈在于：
- 密文多项式的系数需要用大整数表示（数千位）
- 基转换（CRT 重建 → 分解）需要 O(N·log²Q) 的大整数运算
- 每次乘法后都需要完整的 CRT 重建

本方案通过 RNS 表示完全避免了大整数运算，所有中间值都是固定位宽的 uint64_t。

### 2.2 核心数学框架

#### 2.2.1 多项式环与 RNS 表示 (论文 §2.1-2.2)

- **多项式环**: `R = Z[X]/(X^N + 1)`，其中 `N` 是 2 的幂
- **系数环**: `R_q = R / qR = Z_q[X]/(X^N + 1)`
- **RNS 表示**: 给定互素模数集合 `B = {q_0, ..., q_l}`，多项式 `a ∈ R` 的 RNS 表示为：
  ```
  [a]_B = (a mod q_0, a mod q_1, ..., a mod q_l)
  ```
  其中每个 `a mod q_j ∈ R_{q_j}`。

#### 2.2.2 CKKS 编码 (论文 §2.1)

- **Canonical Embedding**: `σ: R → C^{N/2}`，`a(X) ↦ (a(ζ^{e_0}), ..., a(ζ^{e_{N/2-1}}))`
  - `ζ = exp(2πi / 2N)` 是 2N 次本原单位根
  - `e_j = 5^j mod 2N` 是 `Z_{2N}^*` 的生成元的幂次
- **编码**: `Encode(z) = ⌊Δ · σ^{-1}(z)⌉ ∈ R`
  - `Δ` 是缩放因子 (论文中通常 `Δ = 2^p`)
- **解码**: `Decode(m) = σ(m) / Δ`

#### 2.2.3 同态运算 (论文 §2.1, §4)

| 运算 | 公式 | 密文变化 |
|------|------|----------|
| 加法 Add | `ct_add = ct_1 + ct_2` | level 不变 |
| 乘法 Mult | `ct_mult = ct_1 ⊗ ct_2` (含重线性化) | level 不变，scale ≈ Δ² |
| 重缩放 RS | `ct_rs = q_l^{-1} · (ct - ct mod q_l)` | level 减 1 |

#### 2.2.4 Fast Basis Conversion (论文 §2.3)

给定 `a ∈ R` 的 `[a]_B` 表示，转换到基 `C`：
```
Conv_{B→C}([a]_B)_k = [Σ_j (a_j · q̂_j^{-1} mod q_j) · q̂_j] mod p_k
```
其中 `q̂_j = Π_{i≠j} q_i`。

**关键性质 (Lemma 1)**: 输出 `a'` 满足 `a' = a + Q·v`，其中 `|v| < |B|/2`。这个"近似"在后续 ModUp/ModDown 中被控制。

#### 2.2.5 ModUp / ModDown (论文 §3.2)

**Algorithm 1 — ModUp_{C→D}** (论文 §3.2):
```
输入: [a]_C (C = {q_0,...,q_l})
输出: [a']_D (D = B ∪ C)

步骤:
  1. a'_B = Conv_{C→B}([a]_C)     // 用 fast_convert 得到 B 分量
  2. a' = (a'_B, a_C)               // 拼接 B 和 C 分量
```

**Algorithm 2 — ModDown_{D→C}** (论文 §3.2):
```
输入: [b]_D (D = B ∪ C)
输出: [b']_C ≈ [P^{-1} · b]_C

步骤:
  1. a_C = Conv_{B→C}([b]_B)       // 从 B 分量恢复近似的 a
  2. b'_j = P^{-1} · (b_j - a_j) mod q_j  // 减 a 后除以 P
```

#### 2.2.6 Key Switching (论文 §4 KSGen)

评估密钥 `evk = (b, a) ∈ R_{P·Q}^2` 满足：
```
b + a·s ≈ P·s' + e (mod P·Q)
```
其中 `s'` 是目标秘密（如 `s²`, `s̄`, `s_rot`）。

使用方式：对于密文分量 `d_2`（需要 `s² → s` 的转换），计算：
```
(d_0', d_1') = (d_2 · evk.b, d_2 · evk.a)  // 在 D = B∪C 基上
→ ModDown_{D→C} → (d_0'', d_1'')            // 近似除以 P
最终: d_0'' + d_1''·s ≈ d_2·s²              // 重线性化完成
```

#### 2.2.7 Approximate Basis (论文 §3.1)

论文允许密文模数 `q_j` 不完全相等，只需 `|q_j - q| < ε`（`ε` 较小）。这简化了参数选择，但引入了额外的 rescale 误差：
```
|RS误差| ≈ |Δ·m · (Δ/q_l - 1)|
```
当所有 `q_j ≈ Δ` 时，此项可忽略。

### 2.3 噪声分析概要 (论文附录 A)

每次同态操作后噪声增长的近似界：

| 操作 | 噪声增长 |
|------|----------|
| 加密 | `V_enc ≈ Δ·B_m + B_err` |
| 加法 | `V_add ≈ V_1 + V_2` |
| 乘法 | `V_mult ≈ Δ·V_1·V_2 + B_key` |
| Rescale | `V_rs ≈ V/q_l + N/2` |
| ModUp | 额外误差 ≈ `|B|/2` |
| ModDown | 额外误差 ≈ `|B|/2` |

---

## 3. 项目结构详解

### 3.1 目录树

```
CKKSRNS_Demo/
├── README.md               # 项目说明、编译运行指引
├── 2018-931.pdf            # 论文原文 (SAC 2018)
├── demo.py                 # Python 教学演示 (本仓库的教学补充)
│
├── src/                    # C++ 源代码
│   ├── Common.h            # [头文件] 公共标准库引用
│   │
│   ├── Context.h           # [核心] 全局上下文头文件
│   ├── Context.cpp         # [核心] 全局上下文实现 (~1555行，最大文件)
│   │
│   ├── Scheme.h            # [核心] 方案层头文件
│   ├── Scheme.cpp          # [核心] 方案层实现 (加密/解密/同态运算)
│   │
│   ├── SecretKey.h         # 秘密密钥类
│   ├── SecretKey.cpp       # 秘密密钥实现 (HWT 采样)
│   │
│   ├── Ciphertext.h        # 密文数据结构
│   ├── Ciphertext.cpp      # 密文构造与拷贝
│   │
│   ├── Plaintext.h         # 明文数据结构
│   ├── Plaintext.cpp       # 明文构造与拷贝
│   │
│   ├── Key.h               # 评估密钥数据结构
│   ├── Key.cpp             # 评估密钥构造
│   │
│   ├── SchemeAlgo.h        # 高阶算法头文件
│   ├── SchemeAlgo.cpp      # 高阶算法实现 (幂/逆/指数/sigmoid)
│   │
│   ├── TestScheme.h        # 测试入口头文件
│   ├── TestScheme.cpp      # 测试入口实现 (~14种测试用例)
│   │
│   ├── Numb.h              # 数论工具函数头文件
│   ├── Numb.cpp            # 模运算/Barrett/素数检测/单位根
│   │
│   ├── EvaluatorUtils.h    # 辅助工具头文件
│   ├── EvaluatorUtils.cpp  # 随机数生成/数组旋转
│   │
│   ├── StringUtils.h       # 字符串工具头文件
│   ├── StringUtils.cpp     # 打印数组/对比结果
│   │
│   ├── TimeUtils.h         # 计时工具头文件
│   ├── TimeUtils.cpp       # 计时实现
│   │
│   └── main.cpp            # 程序入口 (调用 TestScheme 的各种测试)
│
├── lib/                    # 编译后的库文件目录
└── run/                    # 可执行文件目录
```

### 3.2 编译依赖

| 依赖 | 说明 |
|------|------|
| GCC (≥4.8) | 必须，因使用 `unsigned __int128` |
| `sys/time.h` | POSIX 计时 (Linux/Mac) |
| 标准 C++ 库 | `<iostream>`, `<complex>`, `<vector>`, `<map>` 等 |

### 3.3 核心调用链路

```
main.cpp
  └─→ TestScheme::test*(logN, L, logp, ...)
        ├─→ Context(logN, logp, L, K)          // 初始化全局参数
        ├─→ SecretKey(context)                  // 生成秘密密钥
        ├─→ Scheme(secretKey, context)          // 初始化方案（自动 KeyGen）
        │     ├─→ addEncKey(sk)                 // 生成公钥
        │     └─→ addMultKey(sk)                // 生成乘法密钥
        ├─→ scheme.encrypt(vals, slots, L)      // 编码 + 加密
        │     ├─→ context.encode(m, vals, ...)  // 论文 §2.1 Encode
        │     └─→ encryptMsg(msg)               // 论文 §4 Enc
        ├─→ scheme.mult(ct1, ct2)               // 论文 §4 Mult
        │     ├─→ NTT 域逐点乘 (Karatsuba)
        │     ├─→ raiseAndEqual()               // 论文 §3.2 ModUp
        │     ├─→ mulKey()                       // 乘评估密钥
        │     └─→ backAndEqual()                // 论文 §3.2 ModDown
        ├─→ scheme.reScaleByAndEqual(ct, dl)    // 论文 §4 RS
        └─→ scheme.decrypt(sk, ct)              // 论文 §4 Dec
              ├─→ decryptMsg(sk, ct)
              └─→ context.decode(mx, vals, ...) // 论文 §2.1 Decode
```

---

## 4. 关键代码文件详解

### 4.1 Common.h — 公共头文件

**功能**: 集中引入标准库头文件，所有其他文件通过 `#include "Common.h"` 获得基础类型。

**引入的库**:
- `<iostream>`: 输入输出
- `<stdint.h>`: `uint64_t`, `uint32_t` 等定宽整数
- `<vector>`, `<set>`: 容器
- `<math.h>`: 数学函数
- `<complex>`: 复数支持
- `<sys/time.h>`: 计时

**论文对应**: 无直接对应，是工程基础设施。

---

### 4.2 Numb.h / Numb.cpp — 数论工具库

**功能**: 提供底层模运算原语，是整个方案的数学基础。

#### 核心函数一览

| 函数 | 功能 | 论文对应 |
|------|------|----------|
| `addMod(r, a, b, m)` | `r = (a+b) mod m` | 基础运算 |
| `subMod(r, a, b, m)` | `r = (a-b) mod m` | 基础运算 |
| `mulMod(r, a, b, m)` | `r = (a·b) mod m`（用 `__int128`） | 基础运算 |
| `mulModBarrett(r, a, b, p, pr, twok)` | Barrett reduction 加速的模乘 | 性能优化 |
| `modBarrett(r, a, m, mr, twok)` | Barrett reduction 模运算 | 性能优化 |
| `invMod(x, m)` | `x^{-1} mod m`（Fermat 小定理） | §2.2 基转换需要 |
| `powMod(x, y, p)` | `x^y mod p`（快速幂） | §4 KeyGen |
| `bitReverse(x)` | 32 位整数位反转 | NTT 蝶形排序 |
| `findPrimitiveRoot(m)` | 找模 m 的原根 | NTT 需要 M 次单位根 |
| `findMthRootOfUnity(M, p)` | 找模 p 的 M 次单位根 | NTT 需要 |
| `primeTest(p)` | Miller-Rabin 素数检测 | §4 Setup 生成模数 |
| `gcd(a, b)` | 最大公约数 | 辅助 |

#### Barrett Reduction 原理 (论文性能优化的关键)

Barrett reduction 用预计算的 `μ = ⌊2^{2k} / p⌋` 替代除法：
```
q = ⌊a · μ / 2^{2k}⌋    // 用移位代替除法
r = a - q · p              // 减法
if r ≥ p: r -= p          // 修正
```
C++ 代码中 `pr` 就是 `μ`，`twok` 就是 `2k`。

---

### 4.3 Context.h / Context.cpp — 全局上下文（最大最核心文件）

**功能**: 承载所有参数初始化、预计算表、编码/解码、NTT/INTT、采样、RNS 操作。

#### 4.3.1 构造函数 `Context::Context(logN, logp, L, K, h, sigma)` (论文 §4 Setup)

**参数含义**:
| 参数 | 含义 | 论文对应 |
|------|------|----------|
| `logN` | `log₂(N)`，环维度对数 | §2.1: N 是 2 的幂 |
| `logp` | `log₂(Δ)`，缩放因子对数 | §2.1: Δ = 2^p |
| `L` | 最大 level（密文模数链长度） | §4: Setup 参数 |
| `K` | 特殊模数个数 (通常 K = L+1) | §3.2: B 的大小 |
| `h` | 秘密密钥汉明重量 | §4 KeyGen: χ_key 分布 |
| `sigma` | 误差分布标准差 | §4: χ_err 分布，默认 3.2 |

**构造函数的核心步骤**:

1. **生成密文模数 qVec[0..L-1]** (论文 §3.1 Approximate Basis):
   - `q_0` 约 `2^61` (特殊，稍大)
   - `q_1, ..., q_{L-1}` 约 `2^logp`
   - 所有 `q_j ≡ 1 (mod 2N)` 以支持 NTT
   - 使用 `primeTest()` 逐一验证

2. **生成特殊模数 pVec[0..K-1]** (论文 §3.2):
   - 所有 `p_i ≡ 1 (mod 2N)`
   - 约 `2^logp` 大小
   - 搜索从密文模数之后继续，保证互素

3. **预计算 NTT 相关表**:
   - `qRoots[i]`: 模 `q_i` 的 M 次本原单位根
   - `qRootPows[i][j]`: 单位根的幂次表（蝶形运算用）
   - `qRootScalePows[i][j]`: 缩放后的幂次（Montgomery 形式）
   - `NInvModq[i]`: `N^{-1} mod q_i`（逆 NTT 归一化）
   - `qHatModq[l][i]`, `qHatInvModq[l][i]`: 基转换的预计算系数

4. **预计算模逆/模乘表**:
   - `PModq[i]`, `PInvModq[i]`: `P mod q_i` 及其逆（用于 evalAndEqual）
   - `QModp[i][k]`, `QInvModp[i][k]`: `Q_i mod p_k` 及其逆
   - `qHatModp[l][i][k]`, `pHatModq[k][i]`: 基转换交叉表
   - `qInvModq[i][j]`: `q_i^{-1} mod q_j`（用于 rescale）

5. **预计算旋转群**:
   - `rotGroup[i] = 5^i mod 2N`：`Z_{2N}^*` 的循环子群
   - `ksiPows[j] = exp(2πi·j/2N)`：复数单位根表（用于编码/解码的 FFT）

6. **预计算常数多项式**:
   - `p2coeff`: `Δ² mod q_i`（用于 `multByConst`）
   - `pccoeff`: 特定常数（用于 sigmoid）
   - `p2hcoeff`: `Δ·(Δ/2) mod q_i`

#### 4.3.2 编码/解码 (论文 §2.1)

**`Context::encode(uint64_t* a, complex<double>* v, long slots, long l)`**:
```
1. fftSpecialInv(v, slots)         // 论文 §2.1: σ^{-1} 的 FFT 实现
2. v[j] *= p                        // 乘缩放因子 Δ
3. 将实部/虚部填入多项式系数:
   - mi[idx] = real(v[j])           // 实部放在位置 idx
   - mi[jdx] = imag(v[j])           // 虚部放在位置 N/2 + idx
4. qiNTTAndEqual(mi, i)             // 转 NTT 表示
5. 重复 L 次（每个模数一份）
```

**`Context::decode(uint64_t* a, complex<double>* v, long slots, long l)`**:
```
1. qiINTTAndEqual(tmp, 0)           // 从 NTT 转回系数
2. 提取系数:
   - mir = tmp[idx] / p             // 除以缩放因子
   - mii = tmp[jdx] / p
3. fftSpecial(v, slots)             // 论文 §2.1: σ 的正变换
```

#### 4.3.3 NTT / INTT (论文 §2.2 实现)

**`qiNTTAndEqual(a, index)`**: Cooley-Tukey 蝶形 NTT
```
for m = 1, 2, 4, ..., N/2:        // 蝶形层数 = log₂N
  for i = 0..m-1:                   // 每组蝶形
    W = qRootScalePows[index][m+i]  // 旋转因子 (twiddle factor)
    for j in group:                  // 蝶形运算
      T = a[j+t]
      U = T * W                      // Montgomery 乘法
      a[j+t] = a[j] - U
      a[j] = a[j] + U
```

**`qiINTTAndEqual(a, index)`**: Gentleman-Sande 逆 NTT
```
for m = N, N/2, ..., 2:           // 反向蝶形
  for i = 0..m/2-1:
    W = qRootScalePowsInv[index][m/2+i]
    T = a[j] - a[j+t]
    a[j] += a[j+t]
    a[j+t] = T * W                 // 乘以逆旋转因子
最后: a[j] *= N^{-1} mod q         // 归一化
```

#### 4.3.4 基转换与 ModUp/ModDown (论文 §2.3, §3.2)

**`Context::raiseAndEqual(a, l)`** — 对应论文 §3.2 Algorithm 1 ModUp:
```
1. INTTAndEqual(a, l)               // 转回系数表示
2. 对每个系数乘以 qHatInvModq       // 准备基转换
3. 计算 B 分量:
   rak[n] = Σ_i tmp_i[n] * qHatModp[i][k]  // 论文公式
4. NTT B 分量
```

**`Context::backAndEqual(a, l)`** — 对应论文 §3.2 Algorithm 2 ModDown:
```
1. INTTAndEqual(a, l, K)            // 全部转回系数
2. 对 B 分量乘 pHatInvModp          // 准备基转换
3. 计算 C 分量上的近似值:
   res[n] = Σ_k tmp_k[n] * pHatModq[k][i]  // 从 B 转到 C
4. 减原值并除以 P:
   res = (a_C - res) * P^{-1} mod q_i
5. NTT 结果
```

#### 4.3.5 Rescale (论文 §4 RS)

**`Context::reScaleAndEqual(a, l)`**:
```
1. al = a + (l-1)*N                 // 取最高层模数上的分量
2. qiINTTAndEqual(al, l-1)           // 转回系数表示
3. for i = 0..l-2:                   // 对每个剩余模数
   rai = al mod q_i                  // 取 al 在 q_i 上的值
   qiNTTAndEqual(rai, i)             // 转 NTT
   for n:
     rai[n] = (a_i[n] - rai[n]) * q_l^{-1} mod q_i
   // 等价于: (a - a mod q_l) / q_l
4. 丢弃最高层，l 减 1
```

#### 4.3.6 旋转与共轭 (论文 §4 LeftRotate / Conjugate)

**`Context::leftRot(res, a, l, rotSlots)`**:
```
1. INTT(a)                           // 转回系数
2. for n:
     shift = n * rotGroup[rotSlots] mod M
     if shift < N: res[shift] = a[n]
     else: res[shift-N] = q_i - a[n]  // X^N = -1
3. NTT(res)
```

**`Context::conjugate(res, a, l)`**:
```
for n: res[n] = a[N-1-n]
// 等价于 X → X^{-1} = X^{2N-1}
```

#### 4.3.7 evalAndEqual (论文 §4 KSGen 中的 P 因子注入)

**`Context::evalAndEqual(a, l)`**:
```
1. INTT(a, l)
2. for i = 0..l-1:
     a_i *= PModq[i]                // 每个分量乘 P mod q_i
3. NTT(a, l)
```
这使得 `a` 的 Q 分量被乘以 `P`，等价于将 `a` 嵌入到 `P·Q` 的大环中。

---

### 4.4 SecretKey.h / SecretKey.cpp — 秘密密钥

**功能**: 秘密密钥 `sk = s`，用小汉明重量多项式表示。

**构造 (论文 §4 KeyGen)**:
```cpp
SecretKey::SecretKey(Context& context) {
    sx = new uint64_t[N * (L + K)]();   // 分配 D = B∪C 大小的空间
    context.sampleHWT(sx, L, K);         // 论文 χ_key: 汉明重量 h 的稀疏三元多项式
    context.NTTAndEqual(sx, L, K);       // 全部转 NTT 表示
}
```

**论文对应**: §4 KeyGen，`s ← χ_key`，其中 `χ_key` 是汉明重量为 `h` 的稀疏分布。

---

### 4.5 Key.h / Key.cpp — 评估密钥

**功能**: 评估密钥 `(ax, bx)` 的二元组存储。用于 key switching。

**论文对应**: §4 KSGen 的输出。每个评估密钥满足 `bx + ax·s ≈ P·s' + e (mod P·Q)`。

---

### 4.6 Ciphertext.h / Ciphertext.cpp — 密文

**数据结构**:
```cpp
class Ciphertext {
    uint64_t* bx;     // 密文分量 c0 (NTT 表示)
    uint64_t* ax;     // 密文分量 c1 (NTT 表示)
    long N;           // 环维度
    long slots;       // 打包的 slot 数
    long l;           // 当前 level (使用 l 个密文模数)
};
```

**内存布局**: `bx` 和 `ax` 各占 `N * l` 个 uint64_t，按模数顺序排列：
```
bx = [bx mod q_0 | bx mod q_1 | ... | bx mod q_{l-1}]
ax = [ax mod q_0 | ax mod q_1 | ... | ax mod q_{l-1}]
```

**论文对应**: §4 密文 `ct = (c_0, c_1) ∈ R_{Q_l}^2`。

**注意**: C++ 代码中 `ax` 对应 `c_1`，`bx` 对应 `c_0`（与一些文献的命名不同）。解密时计算 `bx + ax·s`。

---

### 4.7 Plaintext.h / Plaintext.cpp — 明文

**数据结构**:
```cpp
class Plaintext {
    uint64_t* mx;     // 明文多项式 (NTT 表示)
    long N;           // 环维度
    long slots;       // slot 数
    long l;           // level
};
```

**论文对应**: §2.1 Encode 的输出 `m ∈ R`。

---

### 4.8 Scheme.h / Scheme.cpp — 方案核心（第二核心文件）

**功能**: 承载论文 §4 的所有核心算法。

#### 4.8.1 密钥生成 (论文 §4 KeyGen)

**`Scheme::addEncKey(sk)`** — 论文 §4 公钥生成:
```
1. a ← R_Q                             // 均匀随机
2. e ← χ_err                           // 高斯误差
3. b = e - a·s                         // pk = (b, a)
```
存储到 `keyMap[ENCRYPTION]`。

**`Scheme::addMultKey(sk)`** — 论文 §4 KSGen(s², s):
```
1. sxsx = s²                            // s 的平方
2. evalAndEqual(sxsx)                   // 乘 P (论文 KSGen 的 P 因子)
3. e ← χ_err
4. ex = sxsx + e
5. a ← R_{P·Q}                         // 均匀随机
6. b = ex - a·s                         // evk = (b, a)
```
存储到 `keyMap[MULTIPLICATION]`。

**`Scheme::addConjKey(sk)`** — 论文 §4 KSGen(s̄, s):
```
类似 addMultKey，但 s' = conjugate(s)
```

**`Scheme::addLeftRotKey(sk, rot)`** — 论文 §4 KSGen(s_rot, s):
```
类似 addMultKey，但 s' = rotate(s, rot)
```

#### 4.8.2 加密/解密 (论文 §4 Enc / Dec)

**`Scheme::encryptMsg(Plaintext)`** — 公钥加密:
```
1. v ← χ_enc                            // {0, ±1} 采样
2. e0 ← χ_err
3. ax = v·pk.a + e0
4. bx = v·pk.b + e1 + m
// ct = (bx, ax)，解密: bx + ax·s ≈ m
```

**`Scheme::decryptMsg(sk, ct)`**:
```
mx = ax·s + bx                          // 论文 §4 Dec
```
注意 `ax` 对应 `c_1`，`bx` 对应 `c_0`。

#### 4.8.3 同态加法 (论文 §4 Add)

**`Scheme::add(ct1, ct2)`**:
```
res.ax = ct1.ax + ct2.ax                // 逐 NTT 分量相加
res.bx = ct1.bx + ct2.bx
```

#### 4.8.4 同态乘法 (论文 §4 Mult) — 最复杂操作

**`Scheme::mult(ct1, ct2)`** — 使用 Karatsuba 技巧:
```
设 ct1 = (a1, b1), ct2 = (a2, b2)

// 三项乘积 (论文 §4 Mult 的核心)
u0 = a1·a2                               // 对应 s² 项
u1 = a1·b2 + b1·a2                       // 对应 s 项 (Karatsuba 优化)
u2 = b1·b2                               // 对应 1 项

// Karatsuba: 用 (a1+b1)(a2+b2) - u0 - u2 代替 u1

// 重线性化 (论文 §4 Mult + KSGen)
raise(u0)                                // ModUp: u0 从 Q 扩展到 P·Q
(u0_a, u0_b) = (u0·evk.ax, u0·evk.bx)   // 乘评估密钥
back(u0_a); back(u0_b)                   // ModDown: 缩回 Q
res.ax = u0_a + (a1+b1)(a2+b2) - u2 - u0
res.bx = u0_b + u2
```

#### 4.8.5 重缩放 (论文 §4 RS)

**`Scheme::reScaleByAndEqual(ct, dl)`**:
```
for i = 0..dl-1:
    context.reScaleAndEqual(ct.ax, ct.l)  // 论文 Algorithm: RS
    context.reScaleAndEqual(ct.bx, ct.l)
    ct.l -= 1                             // level 减 1
```

#### 4.8.6 旋转 (论文 §4 LeftRotate)

**`Scheme::leftRotateFast(ct, rotSlots)`**:
```
1. bxrot = leftRot(bx, rotSlots)          // 直接旋转 bx
2. bx_rot = leftRot(ax, rotSlots)         // 旋转 ax
3. raise(bx_rot)                          // ModUp
4. (res.ax, res.bx) = bx_rot * evk        // key switching
5. back(res.ax); back(res.bx)             // ModDown
6. res.bx += bxrot                        // 加回旋转后的 bx
```

#### 4.8.7 共轭 (论文 §4 Conjugate)

**`Scheme::conjugate(ct)`**: 类似 leftRotateFast，但用共轭替代旋转。

---

### 4.9 SchemeAlgo.h / SchemeAlgo.cpp — 高阶算法

**功能**: 基于基本同态运算构建高阶算法（论文 §5 的扩展内容）。

#### 核心算法

| 方法 | 功能 | 论文对应 |
|------|------|----------|
| `powerOf2(cipher, logDeg)` | `cipher^{2^logDeg}` (反复平方) | §5 多项式求值 |
| `power(cipher, deg)` | `cipher^deg` (递归分解) | §5 |
| `prodOfPo2(ciphers, logDeg)` | `∏ ciphers` (树状乘法) | §5 |
| `inverse(cipher, steps)` | `1/cipher` (截断几何级数) | §5 函数求值 |
| `exponent(cipher, deg)` | `exp(cipher)` (Taylor 展开) | §5 |
| `sigmoid(cipher, deg)` | `σ(cipher)` (Taylor 展开) | §5 |
| `partialSlotsSum(ct, slots)` | 所有 slot 之和 | §5 聚合操作 |

**`inverse` 的数学原理**:
```
1/(1-x) = 1 + x + x² + x³ + ...
截断到 x^{2^steps}:
  第 0 步: 1 + x
  第 i 步: (1 + x^{2^i}) · res_{i-1}
每次需要 1 次乘法 + 1 次 rescale。
```

---

### 4.10 TestScheme.h / TestScheme.cpp — 测试入口

**功能**: 提供 ~14 种测试函数，覆盖方案的所有核心操作。

| 测试函数 | 测试内容 |
|----------|----------|
| `testEncodeSingle` | 单值编码/解码 |
| `testEncodeBatch` | 批量编码/解码 |
| `testBasic` | 加/乘/常数乘 + rescale |
| `testConjugateBatch` | 共轭 |
| `testimultBatch` | 乘虚数单位 i |
| `testRotateByPo2Batch` | 2 的幂次旋转 |
| `testRotateBatch` | 任意旋转 |
| `testSlotsSum` | slot 求和 |
| `testPowerOf2Batch` | 2 的幂次方 |
| `testPowerBatch` | 任意幂次 |
| `testProdOfPo2Batch` | 树状乘积 |
| `testProdBatch` | 一般乘积 |
| `testInverseBatch` | 求逆 |
| `testLogarithmBatch` | 对数 |
| `testExponentBatch` | 指数 |
| `testSigmoidBatch` | sigmoid |

---

### 4.11 EvaluatorUtils — 辅助工具

| 方法 | 功能 |
|------|------|
| `randomReal(bound)` | 随机实数 |
| `randomComplex(bound)` | 随机复数 |
| `randomCircle(anglebound)` | 单位圆上的随机复数 |
| `randomRealArray(size, bound)` | 随机实数数组 |
| `randomComplexArray(size, bound)` | 随机复数数组 |
| `leftRotateAndEqual(vals, size, rotSize)` | 数组左旋 (用于对比验证) |
| `rightRotateAndEqual(vals, size, rotSize)` | 数组右旋 |

---

### 4.12 辅助工具类

- **StringUtils**: 打印数组、对比明文/密文结果
- **TimeUtils**: 用 `gettimeofday` 计时
- **main.cpp**: 调用 `TestScheme` 的各种测试

---

## 5. 代码与论文精准映射关系

### 5.1 全局映射总表

| 论文位置 | 操作 | C++ 实现位置 |
|----------|------|--------------|
| §2.1 CKKS 编码 | Encode | `Context::encode()` (Context.cpp:458-481) |
| §2.1 CKKS 解码 | Decode | `Context::decode()` (Context.cpp:503-519) |
| §2.2 NTT | Number Theoretic Transform | `Context::qiNTTAndEqual()` (Context.cpp:565-610) |
| §2.2 INTT | Inverse NTT | `Context::qiINTTAndEqual()` (Context.cpp:676-717) |
| §2.3 基转换 | Fast Basis Conversion | `Context::raiseAndEqual()` / `backAndEqual()` |
| §3.1 近似基 | Approximate Basis | `Context` 构造函数中的模数生成逻辑 |
| §3.2 Alg.1 模数提升 | ModUp | `Context::raiseAndEqual()` (Context.cpp:1244-1275) |
| §3.2 Alg.2 模数降低 | ModDown | `Context::backAndEqual()` (Context.cpp:1311-1344) |
| §4 Setup | 参数初始化 | `Context::Context()` (Context.cpp:14-344) |
| §4 KeyGen (pk) | 公钥生成 | `Scheme::addEncKey()` (Scheme.cpp:18-34) |
| §4 KSGen(s²,s) | 乘法密钥 | `Scheme::addMultKey()` (Scheme.cpp:36-59) |
| §4 KSGen(s̄,s) | 共轭密钥 | `Scheme::addConjKey()` (Scheme.cpp:61-83) |
| §4 KSGen(s_rot,s) | 旋转密钥 | `Scheme::addLeftRotKey()` (Scheme.cpp:85-107) |
| §4 Enc | 加密 | `Scheme::encryptMsg()` (Scheme.cpp:163-190) |
| §4 Dec | 解密 | `Scheme::decryptMsg()` (Scheme.cpp:192-198) |
| §4 Add | 同态加法 | `Scheme::add()` (Scheme.cpp:244-252) |
| §4 Sub | 同态减法 | `Scheme::sub()` (Scheme.cpp:263-275) |
| §4 Mult | 同态乘法 | `Scheme::mult()` (Scheme.cpp:295-332) |
| §4 RS | 重缩放 | `Scheme::reScaleByAndEqual()` (Scheme.cpp:628-634) |
| §4 LeftRotate | 左旋转 | `Scheme::leftRotateFast()` (Scheme.cpp:668-690) |
| §4 Conjugate | 共轭 | `Scheme::conjugate()` (Scheme.cpp:765-800) |
| §4 MultByConst | 常数乘 | `Scheme::multByConst()` (Scheme.cpp:524-538) |
| §4 AddConst | 常数加 | `Scheme::addConst()` (Scheme.cpp:476-490) |
| 附录 A 噪声 | 误差采样 | `Context::sampleGauss()` (Context.cpp:1483-1505) |
| §4 χ_key | 秘密密钥采样 | `Context::sampleHWT()` (Context.cpp:1536-1553) |
| §4 χ_enc | 加密随机性 | `Context::sampleZO()` (Context.cpp:1507-1519) |
| 模数生成 | 素数检测 | `primeTest()` (Numb.cpp:173-193) |
| NTT 单位根 | M 次单位根 | `findMthRootOfUnity()` (Numb.cpp:159-170) |
| Barrett 约减 | 快速模运算 | `modBarrett()` (Numb.cpp:35-63) |

### 5.2 加密流程详细映射

```
论文 §4 Enc_pk(m):

  v ← χ_enc                          → Scheme::encryptMsg() 中:
                                        context.sampleZO(vx, Nh, message.l)
                                        context.NTTAndEqual(vx, message.l)

  ct = v·pk + (m+e0, e1)            → 代码:
    ax = v·pk.ax + e0                  context.mul(ax, vx, key.ax, message.l)
    bx = v·pk.bx + e1 + m             context.mul(bx, vx, key.bx, message.l)
                                       context.addAndEqual(bx, message.mx, ...)

  其中 m = Encode(z, Δ)             → Scheme::encode() 中:
                                        context.encode(m, v, slots, l)
```

### 5.3 乘法流程详细映射

```
论文 §4 Mult_{evk}(ct_1, ct_2):

  Step 1: 三项乘积                    → Scheme::mult():
    u0 = a1·a2                          context.mul(axax, cipher1.ax, cipher2.ax)
    u2 = b1·b2                          context.mul(bxbx, cipher1.bx, cipher2.bx)
    u1 = (a1+b1)(a2+b2)-u0-u2          context.add(axbx1, ...); mul; sub...
                                        // Karatsuba 优化

  Step 2: ModUp(u0)                   → context.raiseAndEqual(axax, cipher1.l)
    // axax 从 Q 扩展到 P·Q

  Step 3: 乘评估密钥                   → context.mulKey(axmult, axax, key.ax, ...)
    (d0_D, d1_D) = u0_D · evk         context.mulKey(bxmult, axax, key.bx, ...)

  Step 4: ModDown                      → context.backAndEqual(axmult, ...)
    (d0_C, d1_C) = ModDown(d0_D)      context.backAndEqual(bxmult, ...)

  Step 5: 组合                         → 代码:
    res.ax = d0_C + u1                  context.addAndEqual(axmult, axbx1, ...)
    res.bx = d1_C + u2                  context.addAndEqual(bxmult, bxbx, ...)
```

### 5.4 Rescale 流程详细映射

```
论文 §4 RS_{l,l-1}(ct):

  输入: ct = (c0, c1) with moduli {q_0,...,q_l}

  对 c_i ∈ {c0, c1}:
    1. 取 c^{(l)} 的系数表示            → qiINTTAndEqual(al, l-1)
    2. 对每个 j < l:
       c'^{(j)} = q_l^{-1}·(c^{(j)} - c^{(l)} mod q_j) mod q_j
                                        → for each q_i:
                                            rai = al mod q_i
                                            rai = (ai - rai) * qInvModq[l-1][i]
    3. 丢弃 q_l                         → 分配新的 (l-1)*N 大小的数组

  level 减 1                            → cipher.l -= 1
```

### 5.5 Key Switching 流程映射

```
论文 §4 KSGen(s', s):

  目标: evk = (b, a) 满足 b + a·s = P·s' + e (mod P·Q)

  代码 (以 addMultKey 为例):
    sxsx = s·s                           context.mul(sxsx, sk.sx, sk.sx, L)
    evalAndEqual(sxsx)                   // 乘 P: sxsx → P·s² (mod P·Q)
    e ← χ_err                            context.sampleGauss(ex, L, K)
    ex = sxsx + e
    a ← R_{P·Q}                          context.sampleUniform(ax, L, K)
    b = ex - a·s                          context.mul(bx, ax, sk.sx, L, K)
                                          context.sub2AndEqual(ex, bx, L, K)
```

---

## 6. 核心算法流程图解

### 6.1 完整加密-计算-解密流程

```
明文 z ∈ C^{N/2}
    │
    ▼
[Encode] ─── σ^{-1} + scale ───→ m ∈ R  (论文 §2.1)
    │
    ▼
[Enc] ─── v·pk + (m+e0, e1) ───→ ct = (c0, c1) ∈ R_{Q_L}^2  (论文 §4)
    │
    ├── [Add] ct_add = ct1 + ct2  (论文 §4, level 不变)
    │       │
    ├── [Mult] ct_mult = relinearize(ct1 ⊗ ct2)  (论文 §4)
    │       │    │
    │       │    ├─ [ModUp] u2 → P·Q 基  (论文 §3.2 Alg.1)
    │       │    ├─ [mulKey] u2 · evk  (论文 §4 KSGen)
    │       │    └─ [ModDown] 缩回 Q 基  (论文 §3.2 Alg.2)
    │       │
    │       └── [RS] ct_rs = q_l^{-1}·(ct - ct mod q_l)  (论文 §4)
    │              level 减 1
    │
    ▼
[Dec] ─── c0 + c1·s ───→ m_noisy  (论文 §4)
    │
    ▼
[Decode] ─── σ / scale ───→ z' ∈ C^{N/2}  (论文 §2.1)
    │
    ▼
z' ≈ z (带误差)
```

### 6.2 RNS 模数链与 Level

```
Level L:   {q_0, q_1, ..., q_L}     ← 初始密文 (scale ≈ Δ)
Level L-1: {q_0, q_1, ..., q_{L-1}} ← 1 次 Rescale 后 (scale ≈ Δ)
...
Level 0:   {q_0}                    ← L 次 Rescale 后 (最终解密用 q_0)

特殊模数:  {p_0, ..., p_{K-1}}      ← 仅用于 key switching 中间步骤
           P = p_0 · ... · p_{K-1}
```

### 6.3 乘法中的 RNS 基转换流程

```
          C 基 (密文基)                    D = B ∪ C 基 (扩展基)
          ┌──────────────┐                ┌─────────────────────┐
u2 ∈ R_Q: │ u2 mod q_0    │   ModUp        │ u2 mod p_0  (新)    │
          │ u2 mod q_1    │ ──────────→    │ u2 mod p_1  (新)    │
          │ ...           │  (fast_conv)   │ u2 mod q_0          │
          │ u2 mod q_l    │                │ ...                 │
          └──────────────┘                │ u2 mod q_l          │
                                           └──────────┬──────────┘
                                                      │
                                              mulKey × evk
                                                      │
                                           ┌──────────▼──────────┐
                                           │ relin_D ∈ R_{P·Q}   │
                                           └──────────┬──────────┘
                                                      │
                                                 ModDown
                                               (近似除以 P)
                                                      │
                                           ┌──────────▼──────────┐
                                           │ relin_C ∈ R_Q       │
                                           │ ≈ u2·s² / P · P     │
                                           │ = u2·s²             │
                                           └─────────────────────┘
```

---

## 7. 数学基础速查

### 7.1 符号对照表

| 论文符号 | 代码变量 | 含义 |
|----------|----------|------|
| N | `context.N` | 环维度 (2 的幂) |
| M = 2N | `context.M` | 分圆多项式次数 |
| Δ = 2^p | `context.p` | 缩放因子 |
| L | `context.L` | 最大 level |
| K | `context.K` | 特殊模数个数 |
| q_j | `context.qVec[j]` | 第 j 个密文模数 |
| p_i | `context.pVec[i]` | 第 i 个特殊模数 |
| Q_l = Π_{j≤l} q_j | (隐式) | level l 的模数乘积 |
| P = Π_{i<K} p_i | (隐式) | 特殊模数乘积 |
| ζ = e^{πi/N} | `context.ksiPows[]` | 2N 次本原单位根 |
| h | `context.h` | 秘密密钥汉明重量 |
| σ | `context.sigma` | 误差分布标准差 |
| s | `sk.sx` | 秘密密钥多项式 |
| ct = (c0, c1) | `ct.bx, ct.ax` | 密文 |
| evk = (b, a) | `key.bx, key.ax` | 评估密钥 |

### 7.2 关键恒等式

```
X^N ≡ -1 (mod X^N + 1)                    // 环的约化规则
X^{-1} ≡ X^{N-1} · X^{-N} = -X^{N-1}     // 逆元
q_j ≡ 1 (mod 2N)                          // NTT 友好的素数
ζ^M = 1, ζ^N = -1                         // 单位根性质
[a]_B = (a mod q_0, ..., a mod q_l)       // RNS 表示
Conv_{B→C}([a]_B) = a + Q·v, |v| < |B|/2 // 基转换的近似性
```

### 7.3 噪声预算

```
典型参数: N = 2^15, L = 11, logp = 55, K = 12
- 每个 q_j ≈ 2^55
- P ≈ 2^{55×12} = 2^{660}
- Q_L ≈ 2^{55×11} = 2^{605}
- 初始噪声 ≈ 2^55 (编码噪声)
- 每次乘法噪声平方后乘 Δ ≈ 2^55
- 每次 Rescale 噪声除以 q_l ≈ 2^55
- 最终可支持约 L 层乘法深度
```

---

## 8. 关键数据结构

### 8.1 密文在内存中的布局

```
Ciphertext.ct (l = 3 的示例):

bx (c0):
┌──────────┬──────────┬──────────┐
│ mod q_0  │ mod q_1  │ mod q_2  │
│ N 个系数 │ N 个系数 │ N 个系数 │
└──────────┴──────────┴──────────┘
  N 个 uint64   N 个 uint64   N 个 uint64

ax (c1):
┌──────────┬──────────┬──────────┐
│ mod q_0  │ mod q_1  │ mod q_2  │
│ N 个系数 │ N 个系数 │ N 个系数 │
└──────────┴──────────┴──────────┘

总大小: 2 × l × N × 8 字节
```

### 8.2 Context 预计算表大小估算

| 表 | 大小 | 用途 |
|----|------|------|
| `qRootPows[L][N]` | L·N·8 B | NTT 旋转因子 |
| `qHatModq[L][L]` | L²·8 B | 基转换系数 |
| `qHatModp[L][L][K]` | L²·K·8 B | 跨基转换 |
| `ksiPows[2N+1]` | (2N+1)·16 B | 编码 FFT 用复数根 |

当 N=2^16, L=26, K=27 时:
- `qRootPows`: 26×65536×8 ≈ 13.6 MB
- `qHatModp`: 26²×27×8 ≈ 146 MB
- 总计约 200-300 MB

### 8.3 密钥存储结构

```
Scheme 类内部:
  keyMap:
    [ENCRYPTION=0]     → Key(pk.ax, pk.bx)         // 公钥
    [MULTIPLICATION=1] → Key(mult.ax, mult.bx)     // 乘法密钥
    [CONJUGATION=2]    → Key(conj.ax, conj.bx)     // 共轭密钥

  leftRotKeyMap:
    [1]    → Key(rot_1.ax, rot_1.bx)               // 左旋 1
    [2]    → Key(rot_2.ax, rot_2.bx)               // 左旋 2
    [4]    → Key(rot_4.ax, rot_4.bx)               // 左旋 4
    ...
    [N/4]  → Key(rot_{N/4}.ax, rot_{N/4}.bx)       // 左旋 N/4
```

---

## 附录 A: 编译与调试指南

### A.1 GCC 特有特性

```cpp
unsigned __int128    // 128 位无符号整数 (GCC 扩展)
```
用于 NTT 蝶形运算中的中间乘积，避免溢出。

### A.2 Barrett Reduction 参数

```cpp
qTwok[i] = 2 * (log2(qVec[i]) + 1);  // Barrett 位移量
qrVec[i] = (1 << qTwok[i]) / qVec[i]; // Barrett 预计算常数
```

### A.3 Montgomery 形式

NTT 旋转因子以 Montgomery 形式存储 (`qRootScalePows`)，用于避免 NTT 蝶形中的大数除法。

---

> **文档版本**: v1.0  
> **对应论文**: eprint 2018/931 (SAC 2018)  
> **对应代码版本**: 仓库当前 HEAD
