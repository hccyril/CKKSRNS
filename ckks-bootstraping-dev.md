# 经典 CKKS Bootstrapping 开发文档

## 1. 文档目标与工程边界

本文档面向准备在当前 FullRNS-HEAAN 风格代码库中实现经典 CKKS bootstrapping 的软件工程师。目标是在开始写代码前，明确 CKKS bootstrapping 的数学含义、算法步骤、关键公式、模块划分、接口设计、参数规划、误差预算和测试方法。

当前代码库已经提供 CKKS/HEAAN 的核心能力，包括：

- RNS 模数链、NTT/INTT、CRT 基础运算；
- 编码、解码、加密、解密；
- 同态加法、乘法、常数乘、明文向量乘；
- rescale、modDown、key switching；
- slot rotation、conjugation；
- 基础函数求值测试框架。

Bootstrapping 应作为现有 `Scheme` / `Context` / `SchemeAlgo` 体系上的增量扩展，而不是重写底层同态库。实现的核心目标是：给定一个低层级、低剩余模数的 CKKS 密文 `ct`，同态地近似执行一次“解密再重新编码”的流程，使输出密文 `ct'` 在保持明文近似值的同时恢复到较高 level，从而继续支持后续同态计算。

本文所说“经典 CKKS bootstrapping”主要指 CKKS 原始 bootstrapping 思路中的四段流程：

1. Modulus raising；
2. CoeffToSlot；
3. EvalMod；
4. SlotToCoeff。

## 2. 子任务 1：开发文档大纲

### 2.1 背景与符号

需要先解释 CKKS 的环、模数、scale、编码、密文语义、level 和噪声模型。

### 2.2 Bootstrapping 总流程

需要给出输入输出约束，并把流程拆成 `ModRaise -> CoeffToSlot -> EvalMod -> SlotToCoeff -> normalize`。

### 2.3 算法细节

需要分别说明：

- Modulus raising 如何把低 level 密文提升到 bootstrapping 扩展模数链；
- CoeffToSlot 如何同态执行近似解码线性变换；
- EvalMod 如何近似计算模约减函数；
- SlotToCoeff 如何同态执行近似编码线性变换；
- 每一步的 level、scale、旋转密钥和误差消耗。

### 2.4 工程模块

需要定义：

- bootstrapping 参数结构；
- bootstrapping key 生成；
- 线性变换矩阵生成；
- baby-step giant-step 线性变换求值器；
- polynomial approximator；
- bootstrap pipeline；
- 测试和 benchmark。

### 2.5 开发计划

需要列出从最小可用版本到优化版本的里程碑，以及每个阶段的验收标准。

## 3. CKKS 基础回顾

### 3.1 环与密文空间

令多项式环为：

```text
R = Z[X] / (X^N + 1)
```

其中 `N` 是 2 的幂，通常令 `M = 2N`。模数链由若干 NTT-friendly 素数组成：

```text
Q_L = q_0 q_1 ... q_L
```

RNS 实现中，一个多项式 `a in R_Q` 被表示为：

```text
a = (a mod q_0, a mod q_1, ..., a mod q_L)
```

当前代码中的 `Context::qVec`、`pVec`、`NTT`、`INTT`、`raise`、`back`、`reScale`、`modDown` 等成员与函数就是这部分的实现基础。

### 3.2 CKKS 编码语义

CKKS 把复数向量 `z in C^s` 编码为环元素：

```text
m = Encode(z; Delta) in R_Q
```

其中：

- `s` 是 slot 数，一般 `s <= N/2`；
- `Delta = 2^logp` 是 scale；
- 解码满足 `Decode(m / Delta) ~= z`。

加密后的密文可以写成二元组：

```text
ct = (b, a)
```

以当前代码命名看，`Ciphertext` 内部为 `ax` 和 `bx`。解密语义可以抽象为：

```text
b + a * s = m + e mod Q
```

其中 `s` 是 secret key，`e` 是噪声。不同实现可能使用 `(a, b)` 或 `(b, a)` 的符号顺序，工程实现时以 `Scheme::decryptMsg` 中实际组合方式为准。

### 3.3 Level、Scale 与 Rescale

CKKS 乘法后 scale 近似相乘：

```text
scale(ct1 * ct2) ~= scale(ct1) * scale(ct2)
```

若两个输入 scale 均为 `Delta`，乘法后 scale 为 `Delta^2`，需要 rescale：

```text
ct' = round(ct / q_l)
scale(ct') ~= Delta^2 / q_l
```

通常选择 `q_l ~= Delta`，使 rescale 后 scale 回到 `Delta`。Bootstrapping 的目标就是当 level 即将耗尽时，通过同态执行解密相关变换恢复可用 level。

## 4. Bootstrapping 问题定义

### 4.1 输入

输入密文：

```text
ct_low encrypts m = Delta * z + e mod q
```

其中：

- `ct_low.l` 很小，通常只剩 1 到 2 个普通模数；
- `scale(ct_low) ~= Delta`；
- 明文 slot 值 `z` 需要处于 EvalMod 可处理区间；
- 已生成 conjugation key 与必要 rotation keys；
- bootstrapping 上下文中存在足够长的扩展模数链。

### 4.2 输出

输出密文：

```text
ct_high encrypts Delta * z + e' mod Q_boot
```

要求：

- `ct_high.l` 恢复到 bootstrapping 后目标 level；
- slot 值近似保持；
- `e'` 低于后续计算可接受阈值；
- scale 回到全库约定的 `Delta` 或文档化的 bootstrapping 输出 scale。

### 4.3 核心思想

CKKS bootstrapping 利用密文中“带噪声的消息多项式”：

```text
v = b + a * s = m + e + q * I
```

其中 `I` 是某个整数环元素。普通解密会在私钥端把 `v mod q` 还原成消息；bootstrapping 则先把密文提升到大模数 `Q`，然后同态地近似移除 `q * I` 这部分，使结果回到与 `m` 对应的编码空间。

## 5. 经典 CKKS Bootstrapping 总流程

推荐将工程入口定义为：

```cpp
Ciphertext SchemeAlgo::bootstrap(Ciphertext& ct, BootstrapParams& params);
```

总流程：

```text
Input:  ct at low level, scale Delta
Output: refreshed ct at high level, scale Delta

1. ct0 = modRaise(ct)
2. ct1 = coeffToSlot(ct0)
3. ct2 = evalMod(ct1)
4. ct3 = slotToCoeff(ct2)
5. ct4 = normalizeScaleAndLevel(ct3)
6. return ct4
```

在 CKKS 中，`CoeffToSlot` 和 `SlotToCoeff` 是线性变换，主要消耗 rotation key 和 constant plaintext multiplication。`EvalMod` 是非线性近似函数求值，主要消耗乘法深度。

## 6. Step 1：Modulus Raising

### 6.1 数学说明

输入低模数密文在 `R_q` 上：

```text
ct = (b, a) mod q
```

满足：

```text
b + a * s = m + e mod q
```

Modulus raising 将密文系数从 `R_q` 映射到较大模数 `R_Q`：

```text
ct_up = (b, a) mod Q
```

这里不是普通的“增加精度”，而是选择中心提升：

```text
[x]_q in (-q/2, q/2]  ->  [x]_Q
```

提升后：

```text
b + a * s = m + e + q * I mod Q
```

后续 EvalMod 的任务就是去掉或近似去掉 `q * I` 对应的整数周期项。

### 6.2 工程实现要点

当前 `Context` 已经有 `raiseAndEqual` / `raise`，但需要确认它是否符合 bootstrapping 所需的中心提升语义。若只是从普通基到特殊基的 key switching 辅助提升，则不应直接复用为 bootstrapping 的 `modRaise`。

建议新增：

```cpp
void Context::modRaiseBoot(uint64_t*& poly, long fromLevel, long toLevel);
void SchemeAlgo::modRaiseBootAndEqual(Ciphertext& ct, const BootstrapParams& params);
```

实现要求：

- 输入只包含 `q_0 ... q_l`；
- 输出包含 bootstrapping 所需 `q_0 ... q_Lboot`；
- 对每个系数做中心代表元解释；
- 在新增模数下补齐对应 residue；
- 输出保持 NTT 表示约定与现有 `Ciphertext` 一致。

### 6.3 验收测试

对随机小多项式 `a`：

```text
a mod q -> center lift -> mod Q -> mod q
```

应恢复原值。对边界值 `q/2 - 1`、`q/2`、`q - 1` 应单独测试中心提升符号。

## 7. Step 2：CoeffToSlot

### 7.1 数学说明

`CoeffToSlot` 同态执行 CKKS 解码中的线性变换。可以把编码看作：

```text
Encode(z) = pi^{-1}(z) * Delta
```

解码近似为：

```text
Decode(m) = pi(m / Delta)
```

其中 `pi` 是从环多项式系数表示到 canonical embedding slot 表示的线性映射。`CoeffToSlot` 要在密文态计算：

```text
ct_slots = HomLinearTransform(ct_coeffs, U)
```

使其 slot 中承载：

```text
slots(ct_slots) ~= q * I + m + e
```

更准确地说，经典 CKKS bootstrapping 会把系数表示拆成共轭对或实虚部两路，以便在 slot 上执行周期函数近似。常见做法是输出两个密文：

```text
ct_real ~= Re(pi(v))
ct_imag ~= Im(pi(v))
```

或通过复数 slot 直接承载待 EvalMod 的复数值。

### 7.2 对角线线性变换

任意 slot 线性变换 `A` 可用对角线方法表示：

```text
A * z = sum_{k in D} diag_k(A) * rot(z, k)
```

其中：

- `rot(z, k)` 表示 slot 左旋 `k`；
- `diag_k(A)` 是第 `k` 条循环对角线；
- `D` 是非零对角线集合。

同态求值形式：

```text
ct_out = sum_{k in D} PlainMult(Rotate(ct_in, k), Encode(diag_k(A)))
```

为了降低旋转数量，应实现 baby-step giant-step：

```text
k = i + j * b
A(z) = sum_j rot( sum_i diag_{i,j} * rot(z, i), j*b )
```

其中 `b` 是 baby step 大小。工程上需要权衡：

- baby steps 数量；
- giant steps 数量；
- rotation key 数量；
- 临时密文内存。

### 7.3 工程模块

建议新增：

```cpp
struct LinearTransformPlan {
    long logSlots;
    long levelStart;
    long levelEnd;
    long babyStep;
    vector<long> rotations;
    vector<Plaintext> diagonals;
};

class BootstrapLinearTransform {
public:
    LinearTransformPlan buildCoeffToSlotPlan(const BootstrapParams& params);
    LinearTransformPlan buildSlotToCoeffPlan(const BootstrapParams& params);
    Ciphertext evaluate(const Ciphertext& ct, const LinearTransformPlan& plan);
};
```

如果暂时不引入类，也可放在 `SchemeAlgo` 中：

```cpp
Ciphertext SchemeAlgo::linearTransform(Ciphertext& ct, LinearTransformPlan& plan);
Ciphertext SchemeAlgo::coeffToSlot(Ciphertext& ct, BootstrapParams& params);
Ciphertext SchemeAlgo::slotToCoeff(Ciphertext& ct, BootstrapParams& params);
```

### 7.4 矩阵生成

`CoeffToSlot` 的矩阵来自 CKKS 编码矩阵的逆方向。工程实现不要在运行时反复生成大矩阵；应在参数初始化阶段生成并缓存对角线明文。

需要缓存：

- `coeffToSlotDiags`；
- `slotToCoeffDiags`；
- 每条 diagonal 对应的明文 scale；
- 每条 diagonal 所需 level；
- BSGS 旋转集合。

### 7.5 Level 与 Scale

线性变换主要执行 plaintext multiplication：

```text
ct * pt_diag
```

若 `pt_diag` 使用 scale `Delta` 编码，乘后 scale 变为 `Delta^2`，通常需要 rescale 一层回到 `Delta`。因此一组线性变换至少消耗 1 层。若采用分块或多级矩阵，可能消耗更多层。

## 8. Step 3：EvalMod

### 8.1 数学目标

ModRaise 后的明文包含：

```text
v = m + e + q * I
```

希望同态计算：

```text
v mod q ~= m + e
```

CKKS 中不能精确计算离散模函数，因此使用周期函数近似。常见归一化形式为：

```text
x = v / q
```

则：

```text
v mod q = q * (x - round(x))
```

函数：

```text
f(x) = x - round(x)
```

在整数附近不连续，不能直接用低阶多项式稳定逼近。经典 CKKS bootstrapping 利用正弦函数在小区间内近似线性：

```text
sin(2 * pi * x) ~= 2 * pi * (x - round(x))
```

因此：

```text
x - round(x) ~= (1 / (2 * pi)) * sin(2 * pi * x)
```

再乘回 `q`：

```text
v mod q ~= (q / (2 * pi)) * sin(2 * pi * v / q)
```

实际工程中会根据输入范围、scale 和误差预算选择：

- 直接近似 `sin(2*pi*x)`；
- 近似 sawtooth 函数；
- 使用 Chebyshev 多项式；
- 使用 double-angle / angle-doubling 降低初始多项式次数。

### 8.2 输入区间约束

为了使近似稳定，需要保证归一化后的消息落在小区间：

```text
x = v / q = I + epsilon
```

其中：

```text
|epsilon| <= K_eval
```

`K_eval` 由明文幅度、噪声和缩放共同决定。若 slot 值太大，`EvalMod` 会跨越不连续点导致刷新失败。开发时应在 `BootstrapParams` 中明确：

```cpp
double messageBound;
double evalModInputBound;
```

并在 debug 模式下通过解密测试统计实际输入范围。

### 8.3 多项式近似

设要逼近的函数为：

```text
g(x) = (1 / (2*pi)) * sin(2*pi*x)
```

在区间 `[-B, B]` 上生成 Chebyshev 近似：

```text
P_d(x) = sum_{i=0}^{d} c_i T_i(x / B)
```

其中 `T_i` 是 Chebyshev 多项式：

```text
T_0(x) = 1
T_1(x) = x
T_{i+1}(x) = 2xT_i(x) - T_{i-1}(x)
```

工程可先实现 power basis 版本：

```text
P_d(x) = c_0 + c_1 x + ... + c_d x^d
```

再优化为 Paterson-Stockmeyer 或 baby-step giant-step 多项式求值。

### 8.4 Angle Doubling 版本

如果先用低阶多项式近似：

```text
s_0 ~= sin(2*pi*x / 2^r)
```

再通过倍角公式恢复：

```text
s_{i+1} = 2 * s_i * sqrt(1 - s_i^2)
```

但 `sqrt` 不适合同态直接计算。工程上更常用的可实现近似是使用三角倍角多项式：

```text
sin(2t) = 2 sin(t) cos(t)
cos(2t) = 1 - 2 sin^2(t)
```

若只维护 `sin`，可以用：

```text
sin(2t)^2 = 4s^2(1 - s^2)
```

这会引入符号与平方根问题。因此第一版建议采用固定 Chebyshev 多项式近似 `g(x)`，先保证正确性，再考虑更复杂的 EvalMod 优化。

### 8.5 工程接口

建议新增：

```cpp
struct EvalModPoly {
    vector<double> coeffs;
    double inputBound;
    long degree;
    long depth;
    long scaleBits;
};

Ciphertext SchemeAlgo::evalMod(Ciphertext& ct, const EvalModPoly& poly);
Ciphertext SchemeAlgo::evalPoly(Ciphertext& ct, const vector<double>& coeffs);
```

第一阶段可以手工配置低阶奇函数多项式，例如：

```text
sin(2*pi*x)/(2*pi) ~= x - (2*pi)^2*x^3/3!/(2*pi) + (2*pi)^4*x^5/5!/(2*pi) - ...
```

即：

```text
g(x) ~= x - (2*pi)^2 / 6 * x^3 + (2*pi)^4 / 120 * x^5 - ...
```

但 Taylor 多项式只适合很小区间。正式参数建议使用离线 Chebyshev 拟合。

## 9. Step 4：SlotToCoeff

### 9.1 数学说明

`SlotToCoeff` 是 `CoeffToSlot` 的逆向线性变换，同态执行 CKKS 编码矩阵方向，把 EvalMod 后 slot 中的近似值重新放回多项式系数表示：

```text
ct_coeff = HomLinearTransform(ct_slots, U^{-1})
```

目标：

```text
ct_coeff encrypts Encode(z; Delta) + e_boot
```

### 9.2 工程复用

`SlotToCoeff` 与 `CoeffToSlot` 可以共用同一个线性变换执行器，只是 plan 不同：

```cpp
LinearTransformPlan coeffToSlotPlan;
LinearTransformPlan slotToCoeffPlan;
```

注意：

- 两个矩阵的 diagonal 明文不同；
- 旋转集合可能不同；
- scale 与 level 消耗类似；
- 若 `CoeffToSlot` 拆分成实部和虚部两路，`SlotToCoeff` 需要将两路重新组合。

## 10. 参数设计

建议新增：

```cpp
struct BootstrapParams {
    long logN;
    long logSlots;
    long logDelta;
    long inputLevel;
    long outputLevel;
    long bootstrapLevel;
    long coeffToSlotLevelBudget;
    long evalModLevelBudget;
    long slotToCoeffLevelBudget;
    long linearTransformBabyStep;
    long evalModDegree;
    double messageBound;
    double evalModInputBound;
    bool useBSGS;
    bool debugDecrypt;
};
```

参数约束：

```text
bootstrapLevel >= coeffToSlotLevelBudget
                + evalModLevelBudget
                + slotToCoeffLevelBudget
                + safetyMargin
```

典型第一版配置建议：

```text
logN                  = 14 or 15
logSlots              = logN - 1
logDelta              = 40
inputLevel            = 1 or 2
coeffToSlotBudget     = 1
evalModBudget         = 6 to 12
slotToCoeffBudget     = 1
safetyMargin          = 2
```

说明：实际安全参数与可用深度需要结合模数大小、误差和安全估计工具重新确认。本文档给的是工程开发起点，不是生产安全参数。

## 11. Bootstrapping Key 设计

Bootstrapping 至少需要：

- relinearization key；
- conjugation key；
- CoeffToSlot 所需 rotation keys；
- SlotToCoeff 所需 rotation keys；
- BSGS baby/giant step rotation keys。

建议新增：

```cpp
struct BootstrapKeySpec {
    vector<long> rotations;
    bool needConjugation;
    bool needMultiplication;
};

BootstrapKeySpec SchemeAlgo::collectBootstrapKeySpec(const BootstrapParams& params);
void Scheme::addBootstrapKeys(SecretKey& sk, const BootstrapKeySpec& spec);
```

当前代码已有：

```cpp
void Scheme::addMultKey(SecretKey& secretKey);
void Scheme::addConjKey(SecretKey& secretKey);
void Scheme::addLeftRotKey(SecretKey& secretKey, long rot);
void Scheme::addLeftRotKeys(SecretKey& secretKey);
```

因此第一版可直接复用这些函数生成所需旋转。为了避免生成全部 rotation key 带来内存压力，正式实现应只生成 `collectBootstrapKeySpec` 收集出的 rotation 集合。

## 12. 模块划分建议

### 12.1 新增文件

建议新增：

```text
src/BootstrapParams.h
src/BootstrapParams.cpp
src/Bootstrapper.h
src/Bootstrapper.cpp
src/LinearTransform.h
src/LinearTransform.cpp
src/PolyApprox.h
src/PolyApprox.cpp
src/TestBootstrap.h
src/TestBootstrap.cpp
```

若希望保持项目风格更简单，也可以先把入口放入 `SchemeAlgo`：

```text
SchemeAlgo::bootstrap
SchemeAlgo::coeffToSlot
SchemeAlgo::evalMod
SchemeAlgo::slotToCoeff
```

但线性变换和多项式近似建议独立文件，因为这两部分后续会持续优化。

### 12.2 模块职责

`BootstrapParams`：

- 保存 bootstrapping 参数；
- 校验 level budget；
- 生成默认参数；
- 打印参数摘要。

`Bootstrapper`：

- 串联完整 pipeline；
- 管理中间 scale 和 level；
- debug 模式下输出误差统计；
- 对外提供高层 API。

`LinearTransform`：

- 生成 CoeffToSlot / SlotToCoeff diagonal；
- 管理 diagonal 明文缓存；
- 执行 naive diagonal method；
- 执行 BSGS method。

`PolyApprox`：

- 存储 EvalMod 多项式系数；
- 支持 Taylor 初版；
- 支持 Chebyshev 离线系数导入；
- 提供 Horner / Paterson-Stockmeyer 求值。

`TestBootstrap`：

- 单元测试；
- 端到端刷新测试；
- level/scale/误差跟踪；
- benchmark。

## 13. API 设计建议

高层 API：

```cpp
class Bootstrapper {
public:
    Bootstrapper(Scheme& scheme, Context& context, const BootstrapParams& params);

    void prepare();
    BootstrapKeySpec keySpec() const;

    Ciphertext bootstrap(Ciphertext& ct);

private:
    Ciphertext modRaise(Ciphertext& ct);
    Ciphertext coeffToSlot(Ciphertext& ct);
    Ciphertext evalMod(Ciphertext& ct);
    Ciphertext slotToCoeff(Ciphertext& ct);
    void normalize(Ciphertext& ct);
};
```

参数校验 API：

```cpp
bool BootstrapParams::validate(string* error) const;
long BootstrapParams::requiredLevel() const;
vector<long> BootstrapParams::requiredRotations() const;
```

测试入口：

```cpp
void TestScheme::testBootstrap(long logN, long logp, long L, long logSlots);
void TestScheme::testCoeffToSlot(long logN, long logp, long L, long logSlots);
void TestScheme::testEvalMod(long logN, long logp, long L, long logSlots);
```

## 14. 关键实现细节

### 14.1 NTT 表示一致性

现有密文多项式通常在 NTT 域中存储。新增 diagonal 明文、临时多项式、modRaise 输出必须明确：

- 输入是否在 NTT 域；
- 输出是否在 NTT 域；
- multiplication 前是否需要 NTT；
- rotation 前是否要求 NTT。

建议在每个新增函数注释中写清楚 representation contract。

### 14.2 Scale 对齐

每一步结束后记录：

```text
ct.l
ct.slots
approx scale
expected semantic value
```

推荐新增 debug helper：

```cpp
void Bootstrapper::traceState(const string& name, Ciphertext& ct);
```

如果现有 `Ciphertext` 没有显式 scale 字段，开发时要把 scale 作为参数约定并在文档和测试中跟踪。

### 14.3 常数编码

EvalMod 中的常数系数、线性变换 diagonal 都需要编码为 plaintext。要统一：

```text
plaintext scale = Delta
plaintext level = current ct level
```

否则乘法后 scale 和 rescale 会偏离预期。

### 14.4 旋转方向

当前代码有：

```cpp
leftRotateFast
leftRotateByPo2
rightRotateByPo2
leftRotate
rightRotate
```

线性变换 diagonal 公式中的 `rot(z, k)` 必须与实现中的 left/right rotation 对齐。建议先在明文数组上实现相同 diagonal transform，并与密文 transform 对比。

### 14.5 内存管理

当前代码大量使用裸指针 `uint64_t*`。新增模块应遵循现有风格，但建议在新代码内部尽量使用 RAII 包装临时缓冲区，至少做到：

- 每个 `new[]` 都有明确释放路径；
- 中间异常或提前返回不会泄漏；
- BSGS 临时密文数量受参数控制。

如果暂不引入智能指针，测试阶段应配合 AddressSanitizer 或 Valgrind。

## 15. 误差预算

Bootstrapping 总误差可粗略拆为：

```text
e_boot = e_modraise
       + e_coeffToSlot
       + e_evalModApprox
       + e_evalModArithmetic
       + e_slotToCoeff
       + e_rescale
```

其中：

- `e_modraise` 来自中心提升和原始噪声；
- `e_coeffToSlot` 来自线性变换明文编码误差和 rescale；
- `e_evalModApprox` 来自多项式近似误差；
- `e_evalModArithmetic` 来自同态乘法和 relinearization；
- `e_slotToCoeff` 来自逆线性变换；
- `e_rescale` 来自每次 rescale rounding。

端到端测试建议统计：

```text
max_i |z_i - z'_i|
avg_i |z_i - z'_i|
precisionBits = -log2(maxError)
```

验收目标第一版：

```text
maxError < 2^-15
```

优化版可逐步提高到：

```text
maxError < 2^-25 or 2^-30
```

具体目标应根据 `logDelta` 和 EvalMod 多项式阶数调整。

## 16. 测试计划

### 16.1 单元测试

`modRaiseBoot`：

- 小模数人工样例；
- 中心提升边界；
- lift 后再降回原模数一致。

`LinearTransform`：

- 明文矩阵乘法 vs diagonal transform；
- naive diagonal vs BSGS；
- slot rotation 方向；
- CoeffToSlot 后明文解密结果与软件 decode 对齐。

`EvalMod`：

- 明文 double 多项式求值；
- 密文 polynomial evaluation；
- 输入区间边界；
- Taylor 和 Chebyshev 系数误差对比。

### 16.2 集成测试

端到端：

```text
z -> encrypt -> consume levels -> bootstrap -> decrypt -> compare z
```

测试数据：

- 全 0；
- 小常数；
- 随机复数，幅度小于 `messageBound`；
- 接近 EvalMod 边界的输入；
- 实数-only；
- 复数 slot；
- 不同 slot 数。

### 16.3 回归测试指标

每次测试打印：

```text
logN, logSlots, logDelta, L
level before / after
max error
average error
estimated precision bits
runtime per stage
rotation count
```

## 17. 开发里程碑

### 17.1 M1：文档与参数骨架

产出：

- `BootstrapParams`；
- 参数校验；
- bootstrapping key spec 收集；
- 测试入口空壳。

验收：

- 能打印完整参数；
- 能列出所需 rotation keys；
- 不影响现有测试。

### 17.2 M2：线性变换基础版

产出：

- diagonal matrix 生成；
- naive diagonal transform；
- CoeffToSlot / SlotToCoeff 明文对齐测试。

验收：

- 小 `logN` 下明文和密文线性变换一致；
- rotation 方向确认。

### 17.3 M3：EvalMod 基础版

产出：

- Taylor 或手动系数多项式；
- Horner 求值；
- scale/level 跟踪。

验收：

- 小范围输入可同态近似 `x - round(x)`；
- 误差统计稳定。

### 17.4 M4：端到端 Bootstrap MVP

产出：

- `bootstrap(ct)` 串联所有步骤；
- debug trace；
- 端到端测试。

验收：

- 对小幅度随机输入，bootstrap 后可解密回原 slot；
- 输出 level 高于输入；
- precision bits 达到第一版目标。

### 17.5 M5：优化

产出：

- BSGS 线性变换；
- Chebyshev EvalMod；
- rotation key 最小集合；
- benchmark。

验收：

- 相比 naive 版本显著减少旋转和运行时间；
- 精度不低于 MVP；
- 内存占用可解释。

## 18. 风险与注意事项

### 18.1 EvalMod 是最大风险

Bootstrapping 是否成功主要取决于 EvalMod 的输入范围和近似精度。若输入未归一化到多项式可逼近区间，即使其他模块正确也会失败。

### 18.2 当前代码没有显式 scale 字段

如果 `Ciphertext` 只保存 `l` 和 `slots`，则 bootstrapping 实现必须非常严格地按约定维护 scale。建议后续重构时给密文或调试结构增加 scale metadata。

### 18.3 线性变换矩阵容易方向错

CoeffToSlot / SlotToCoeff 的矩阵方向、共轭顺序、bit-reverse 顺序、rotation 正负号都容易出错。必须先用明文模拟器对齐，再进入密文测试。

### 18.4 参数不能只看能跑通

能解密通过不代表安全，也不代表误差预算合理。生产参数需要结合安全估计、模数总 bit 长度、secret distribution、ring dimension 统一评估。

## 19. 推荐代码落点

当前仓库建议按以下方式接入：

```text
src/BootstrapParams.h/.cpp
src/LinearTransform.h/.cpp
src/PolyApprox.h/.cpp
src/Bootstrapper.h/.cpp
src/TestBootstrap.h/.cpp
```

并在现有文件中做少量扩展：

```text
src/Scheme.h/.cpp
  - addBootstrapKeys

src/SchemeAlgo.h/.cpp
  - bootstrap 高层封装，或转调 Bootstrapper

src/Context.h/.cpp
  - modRaiseBoot
  - bootstrapping diagonal/plaintext 生成辅助

src/TestScheme.h/.cpp
  - testBootstrap 入口
```

## 20. 最小伪代码

```text
function Bootstrap(ct):
    assert ct.level <= params.inputLevel
    assert ct.scale ~= Delta

    ct = ModRaiseBoot(ct, params.bootstrapLevel)
    Trace("after mod raise", ct)

    ct = CoeffToSlot(ct, params.coeffToSlotPlan)
    RescaleToDelta(ct)
    Trace("after coeff to slot", ct)

    ct = NormalizeForEvalMod(ct, q)
    ct = EvalPoly(ct, params.evalModPoly)
    ct = ScaleBackAfterEvalMod(ct, q)
    Trace("after eval mod", ct)

    ct = SlotToCoeff(ct, params.slotToCoeffPlan)
    RescaleToDelta(ct)
    Trace("after slot to coeff", ct)

    ct = ModDownTo(ct, params.outputLevel)
    return ct
```

## 21. 第一版实现顺序

建议实际开发时严格按以下顺序推进：

1. 增加 `BootstrapParams` 和 debug trace。
2. 实现明文级 diagonal transform 验证工具。
3. 实现密文 naive linear transform。
4. 跑通 CoeffToSlot / SlotToCoeff 单独测试。
5. 实现 EvalMod 的明文多项式拟合与误差脚本。
6. 实现密文 EvalMod。
7. 串联完整 bootstrapping。
8. 替换 naive linear transform 为 BSGS。
9. 调整参数和多项式阶数，提升精度与速度。

## 22. 完成标准

开发完成后，至少满足：

- `bootstrap(ct)` 对外 API 稳定；
- 所需 rotation keys 可自动收集；
- 端到端测试覆盖实数、复数、随机、小幅度、边界输入；
- 输出密文 level 明确高于输入；
- 输出 scale 与库内约定一致；
- 每个阶段有 runtime 和误差日志；
- 文档中的参数、公式、模块与代码实现保持同步。
