/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_CONTEXT_H_   // 头文件保护宏
#define HEAANNTT_CONTEXT_H_

#include <complex>             // complex<double> 复数类型
#include <chrono>              // 时间相关
#include <map>                 // map 容器（用于 Taylor 系数表）

#include "Common.h"            // 公共头文件
#include "Numb.h"              // 数论工具函数

#define Q0_BIT_SIZE 61         // 第一个密文模数 q_0 的位宽（约 2^61）

using namespace std;

// 预定义函数名称常量（用于 taylorCoeffsMap 的键）
static string LOGARITHM = "Logarithm"; ///< log(x)  对数函数名
static string EXPONENT  = "Exponent"; ///< exp(x)  指数函数名
static string SIGMOID   = "Sigmoid"; ///< sigmoid(x) = exp(x) / (1 + exp(x))  Sigmoid 函数名

// ============================================================================
// 全局上下文类 Context —— 项目中最核心的类
// ============================================================================
// 承载所有参数初始化、预计算表、编码/解码、NTT/INTT、采样、RNS 操作
// 对应论文 §4 Setup: 在构造时完成所有预计算
// 对应论文 §2.1-2.3: 编码/解码/基转换
// ============================================================================
class Context {
public:

	// ========== 加密参数（论文 §4 Setup 输入参数）==========
	long logN; ///< Logarithm of Ring Dimension   log₂(N)，环维度对数
	long logNh; ///< Logarithm of Ring Dimension - 1   log₂(N) - 1 = log₂(N/2)
	long L; ///< Maximum Level that we want to support   最大 level（密文模数链长度）
	long K; ///< The number of special modulus (usually L + 1)   特殊模数个数（通常 K = L+1）

	long N;    // 环维度 N = 2^logN（论文 §2.1: R = Z[X]/(X^N+1)）
	long M;    // M = 2N，分圆多项式 Φ_M(X) = X^N + 1 的次数
	long Nh;   // N/2 = slot 数

	long logp; // log₂(Δ)，缩放因子对数
	long p;    // Δ = 2^logp，缩放因子（论文中记为 Δ）

	long h;        // 秘密密钥的汉明重量（默认 64，论文 §4 χ_key 分布参数）
	double sigma;  // 高斯误差分布的标准差（默认 3.2，论文 §4 χ_err 分布参数）

	// ========== 模数数组 ==========
	uint64_t* qVec;  // 密文模数数组 qVec[0..L-1]，每个 q_j ≡ 1 (mod 2N) 以支持 NTT（论文 §3.1）
	uint64_t* pVec;  // 特殊模数数组 pVec[0..K-1]，每个 p_i ≡ 1 (mod 2N)（论文 §3.2）

	// ========== Barrett Reduction 预计算表 ==========
	uint64_t* qrVec; // Barrett reduction 预计算常数 μ_q[i] = ⌊2^{twok_i}/q_i⌋
	uint64_t* prVec; // Barrett reduction 预计算常数 μ_p[i] = ⌊2^{twok_i}/p_i⌋

	long* qTwok; // Barrett 位移量 twok_q[i] = 2*(⌊log₂(q_i)⌋+1)
	long* pTwok; // Barrett 位移量 twok_p[i] = 2*(⌊log₂(p_i)⌋+1)

	// ========== Montgomery Reduction 预计算表 ==========
	uint64_t* qkVec; // Montgomery 预计算常数（用于 NTT 蝶形中的快速约减）
	uint64_t* pkVec; // Montgomery 预计算常数

	uint64_t* qdVec; // 2*q_i（用于 NTT 中的模约减）
	uint64_t* pdVec; // 2*p_i

	// ========== 模逆预计算表 ==========
	uint64_t* qInvVec; // q_i 的逆元相关预计算值（用于 Montgomery/Barrett reduction 的蝶形运算）
	uint64_t* pInvVec; // p_i 的逆元相关预计算值

	// ========== NTT 单位根 ==========
	uint64_t* qRoots;     // 模 q_i 的 M 次本原单位根 ω_i（论文 §2.2: ω^M ≡ 1 mod q_i）
	uint64_t* pRoots;     // 模 p_i 的 M 次本原单位根

	uint64_t* qRootsInv;  // ω_i 的逆元 ω_i^{-1}（用于 INTT）
	uint64_t* pRootsInv;  // ω_i^{-1}

	// ========== NTT 旋转因子幂次表 ==========
	uint64_t** qRootPows;     // [i][j]: q_i 下 ω_i^j（NTT 蝶形的旋转因子，bit-reversal 排列）
	uint64_t** pRootPows;     // [i][j]: p_i 下 ω_i^j

	uint64_t** qRootPowsInv;  // [i][j]: q_i 下 ω_i^{-j}（INTT 用）
	uint64_t** pRootPowsInv;  // [i][j]: p_i 下 ω_i^{-j}

	// ========== N 的逆元 ==========
	uint64_t* NInvModq;  // [i]: N^{-1} mod q_i（INTT 最后的归一化因子）
	uint64_t* NInvModp;  // [k]: N^{-1} mod p_i

	// ========== NTT 缩放幂次表（Montgomery 形式）==========
	uint64_t** qRootScalePows;      // [i][j]: ω_i^j 的 Montgomery 缩放形式（用于 NTT 蝶形乘法）
	uint64_t** pRootScalePows;      // [i][j]: p_i 版本

	uint64_t** qRootScalePowsOverq; // [i][j]: ⌊(ω_i^j << 64) / q_i⌋（浮点 NTT 优化用）
	uint64_t** pRootScalePowsOverp; // [i][j]: p_i 版本

	uint64_t** qRootScalePowsInv;   // [i][j]: ω_i^{-j} 的 Montgomery 缩放形式（INTT 用）
	uint64_t** pRootScalePowsInv;   // [i][j]: p_i 版本

	uint64_t* NScaleInvModq; // [i]: N^{-1} 的缩放形式（拆成两个 2^32 乘法的 Montgomery 形式）
	uint64_t* NScaleInvModp; // [k]: p_i 版本

	// ========== 基转换预计算表（论文 §2.3 Fast Basis Conversion）==========
	uint64_t** qHatModq;     // [l][i]: q̂_i^(l) mod q_i = (Π_{j≠i, j≤l} q_j) mod q_i
	uint64_t* pHatModp;      // [k]: p̂_k mod p_k = (Π_{j≠k, j<K} p_j) mod p_k

	uint64_t** qHatInvModq;  // [l][i]: (q̂_i^(l))^{-1} mod q_i（基转换中需要乘逆元）
	uint64_t* pHatInvModp;   // [k]: (p̂_k)^{-1} mod p_k

	uint64_t*** qHatModp;    // [l][i][k]: (q̂_i^(l)) mod p_k（C→B 基转换交叉表）

	uint64_t** pHatModq;     // [k][i]: (p̂_k) mod q_i（B→C 基转换交叉表）

	// ========== P 和 Q 的模逆预计算表 ==========
	uint64_t* PModq;         // [i]: P mod q_i = (Π_k p_k) mod q_i（evalAndEqual 中使用）
	uint64_t* PInvModq;      // [i]: P^{-1} mod q_i（ModDown 中使用）

	uint64_t** QModp;        // [i][k]: Q_i mod p_k = (Π_{j≤i} q_j) mod p_k
	uint64_t** QInvModp;     // [i][k]: Q_i^{-1} mod p_k

	uint64_t** qInvModq;     // [i][j]: q_i^{-1} mod q_j（rescale 中需要 q_l^{-1} mod q_j）

	long* rotGroup; ///< precomputed rotation group indexes  旋转群预计算: rotGroup[i] = 5^i mod 2N（论文 §4 左旋转操作使用）

	complex<double>* ksiPows; ///< precomputed ksi powers  复数单位根表: ksiPows[j] = exp(2πi·j/M)（编码/解码的 FFT 使用）

	map<string, double*> taylorCoeffsMap; ///< precomputed taylor coefficients  Taylor 系数表: 按函数名查找多项式系数

	uint64_t* p2coeff;   // [i*N+n]: Δ² mod q_i（multByConst 中用到 Δ²）
	uint64_t* pccoeff;   // 特定常数多项式（sigmoid 近似中使用）
	uint64_t* p2hcoeff;  // [i*N+n]: Δ·(Δ/2) mod q_i（sigmoid 中的 1/2 项）

	// ========== 构造函数（论文 §4 Setup）==========
	// 初始化所有参数、生成模数、预计算所有查找表
	Context(long logN, long logp, long L, long K, long h = 64, double sigma = 3.2);

	// ========== FFT / 编码辅助函数 ==========
	void arrayBitReverse(complex<double>* vals, const long size);  // 复数数组 bit-reversal 置换
	void arrayBitReverse(uint64_t* vals, const long size);         // uint64 数组 bit-reversal 置换

	void fft(complex<double>* vals, const long size);              // 标准 FFT（Cooley-Tukey 蝶形）
	void fftInvLazy(complex<double>* vals, const long size);       // 逆 FFT（不归一化）
	void fftInv(complex<double>* vals, const long size);           // 逆 FFT（含归一化 /size）

	void fftSpecial(complex<double>* vals, const long size);       // CKKS 专用 FFT（使用旋转群 rotGroup 的指数）
	void fftSpecialInvLazy(complex<double>* vals, const long size); // CKKS 专用逆 FFT（不归一化）
	void fftSpecialInv(complex<double>* vals, const long size);    // CKKS 专用逆 FFT（含归一化）

	// ========== 编码 / 解码（论文 §2.1）==========
	void encode(uint64_t* ax, complex<double>* vals, long slots, long l);  // 批量复数编码
	void encode(uint64_t* ax, double* vals, long slots, long l);           // 批量实数编码（未实现）

	void encodeSingle(uint64_t* ax, complex<double>& val, long l);  // 单复数编码
	void encodeSingle(uint64_t* ax, double val, long l);            // 单实数编码（未实现）

	void decode(uint64_t* ax, complex<double>* vals, long slots, long l);  // 批量复数解码
	void decode(uint64_t* ax, double* vals, long slots, long l);           // 批量实数解码（未实现）

	void decodeSingle(uint64_t* ax, complex<double>& val, long l);  // 单复数解码
	void decodeSingle(uint64_t* ax, double val, long l);            // 单实数解码（未实现）

	// ========== NTT / INTT（论文 §2.2）==========
	void qiNTT(uint64_t* res, uint64_t* a, long index);   // 模 q_index 的 NTT（输出到 res）
	void piNTT(uint64_t* res, uint64_t* a, long index);   // 模 p_index 的 NTT

	void NTT(uint64_t* res, uint64_t* a, long l, long k = 0);  // 批量 NTT: 前 l 个用 qVec，后 k 个用 pVec

	void qiNTTAndEqual(uint64_t* a, long index);   // 原地 NTT（模 q_index）
	void piNTTAndEqual(uint64_t* a, long index);   // 原地 NTT（模 p_index）

	void NTTAndEqual(uint64_t* a, long l, long k = 0);  // 批量原地 NTT

	void qiINTT(uint64_t* res, uint64_t* a, long index);  // 模 q_index 的 INTT
	void piINTT(uint64_t* res, uint64_t* a, long index);  // 模 p_index 的 INTT

	void INTT(uint64_t* res, uint64_t* a, long l, long k = 0);  // 批量 INTT

	void qiINTTAndEqual(uint64_t* a, long index);  // 原地 INTT（模 q_index）
	void piINTTAndEqual(uint64_t* a, long index);  // 原地 INTT（模 p_index）

	void INTTAndEqual(uint64_t* a, long l, long k = 0);  // 批量原地 INTT

	// ========== 逐分量环运算（在 RNS 的每个模数上执行）==========
	void qiNegate(uint64_t* res, uint64_t* a, long index);  // 取负（模 q_index）
	void piNegate(uint64_t* res, uint64_t* a, long index);  // 取负（模 p_index）

	void negate(uint64_t* res, uint64_t* a, long l, long k = 0);  // 批量取负

	void qiNegateAndEqual(uint64_t* a, long index);  // 原地取负（模 q_index）
	void piNegateAndEqual(uint64_t* a, long index);  // 原地取负（模 p_index）

	void negateAndEqual(uint64_t* a, long l, long k = 0);  // 批量原地取负

	// 加常数
	void qiAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void piAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void addConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k = 0);
	void qiAddConstAndEqual(uint64_t* a, uint64_t c, long index);
	void piAddConstAndEqual(uint64_t* a, uint64_t c, long index);
	void addConstAndEqual(uint64_t* a, uint64_t c, long l, long k = 0);

	// 减常数
	void qiSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void piSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void subConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k = 0);
	void qiSubConstAndEqual(uint64_t* a, uint64_t c, long index);
	void piSubConstAndEqual(uint64_t* a, uint64_t c, long index);
	void subConstAndEqual(uint64_t* a, uint64_t c, long l, long k = 0);

	// 多项式加法
	void qiAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index);  // a+b mod q_index
	void piAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index);  // a+b mod p_index
	void add(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);  // 批量加法
	void qiAddAndEqual(uint64_t* a, uint64_t* b, long index);  // 原地加法（模 q_index）
	void piAddAndEqual(uint64_t* a, uint64_t* b, long index);
	void addAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);  // 批量原地加法

	// 多项式减法
	void qiSub(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void piSub(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void sub(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);
	void qiSubAndEqual(uint64_t* a, uint64_t* b, long index);
	void piSubAndEqual(uint64_t* a, uint64_t* b, long index);
	void subAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	// sub2: b = a - b（交换减法，用于 Karatsuba 优化中）
	void qiSub2AndEqual(uint64_t* a, uint64_t* b, long index);  // b = a - b mod q_index
	void piSub2AndEqual(uint64_t* a, uint64_t* b, long index);
	void sub2AndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	// 标量乘
	void qiMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index);
	void piMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index);
	void mulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long l, long k = 0);
	void qiMulConstAndEqual(uint64_t* res, uint64_t cnst, long index);
	void piMulConstAndEqual(uint64_t* res, uint64_t cnst, long index);
	void mulConstAndEqual(uint64_t* res, uint64_t cnst, long l, long k = 0);

	// 多项式逐点乘（NTT 域下的环乘 = 逐点乘）
	void qiMul(uint64_t* res, uint64_t* a, uint64_t* b, long index);  // 逐点乘 mod q_index
	void piMul(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void mul(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);  // 批量逐点乘

	void mulKey(uint64_t* res, uint64_t* a, uint64_t* b, long l);  // 与评估密钥的逐点乘（前 l 个用 q，后 K 个用 p）

	void qiMulAndEqual(uint64_t* a, uint64_t* b, long index);  // 原地逐点乘
	void piMulAndEqual(uint64_t* a, uint64_t* b, long index);
	void mulAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	// 平方
	void qiSquare(uint64_t* res, uint64_t* a, long index);
	void piSquare(uint64_t* res, uint64_t* a, long index);
	void square(uint64_t* res, uint64_t* a, long l, long k = 0);
	void qiSquareAndEqual(uint64_t* a, long index);
	void piSquareAndEqual(uint64_t* a, long index);
	void squareAndEqual(uint64_t* a, long l, long k = 0);

	// ========== eval / raise / back（论文 §3.2 / §4 KSGen 的核心操作）==========
	void evalAndEqual(uint64_t* a, long l);   // 乘以 P: a_i *= P mod q_i（论文 KSGen 中的 P 因子注入）

	void raiseAndEqual(uint64_t*& a, long l);  // 论文 §3.2 Algorithm 1 ModUp_{C->D}: 从 C 基扩展到 D = B∪C 基
	void raise(uint64_t* res, uint64_t* a, long l);

	void backAndEqual(uint64_t*& a, long l);   // 论文 §3.2 Algorithm 2 ModDown_{D->C}: 从 D = B∪C 基缩回 C 基（近似除以 P）
	void back(uint64_t* res, uint64_t* a, long l);

	// ========== Rescale（论文 §4 RS）==========
	void reScaleAndEqual(uint64_t*& a, long l);  // 丢弃最高层模数 q_{l-1}，并对剩余模数除以 q_{l-1}
	void reScale(uint64_t* res, uint64_t* a, long l);

	// ========== ModDown（仅丢弃模数，不做除法）==========
	void modDownAndEqual(uint64_t*& a, long l, long dl);  // 丢弃最高的 dl 个模数
	uint64_t* modDown(uint64_t* a, long l, long dl);

	// ========== 旋转与共轭（论文 §4 LeftRotate / Conjugate）==========
	void leftRot(uint64_t* res, uint64_t* a, long l, long rotSlots);       // 左旋转: a(X) → a(X^{5^rotSlots})
	void leftRotAndEqual(uint64_t* a, long l, long rotSlots);

	void conjugate(uint64_t* res, uint64_t* a, long l);                    // 共轭: a(X) → a(X^{-1})
	void conjugateAndEqual(uint64_t* a, long l);

	void mulByMonomial(uint64_t* res, uint64_t* a, long l, long mdeg);     // 乘单项式: a(X) → a(X)·X^{mdeg}
	void mulByMonomialAndEqual(uint64_t* a, long l, long mdeg);

	// ========== 采样函数（论文 §4 中的随机分布）==========
	void sampleGauss(uint64_t* res, long l, long k = 0);  // 离散高斯采样（论文 χ_err 分布，σ=3.2）
	void sampleZO(uint64_t* res, long s, long l, long k = 0);  // {0, ±1} 采样（论文 χ_enc 分布）
	void sampleUniform(uint64_t* res, long l, long k = 0);  // 均匀采样（论文中 a ← R_Q）
	void sampleHWT(uint64_t* res, long l, long k = 0);  // 汉明重量采样（论文 χ_key 分布，h 个非零系数）

};

#endif
