/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAAN_SCHEMEALGO_H_     // 头文件保护宏
#define HEAAN_SCHEMEALGO_H_

#include "Common.h"             // 公共头文件
#include "EvaluatorUtils.h"     // 辅助工具（随机数、数组旋转）
#include "Plaintext.h"          // 明文类
#include "SecretKey.h"          // 秘密密钥类
#include "Ciphertext.h"         // 密文类
#include "Scheme.h"             // 方案核心类

// ============================================================================
// 高阶算法类：基于 Scheme 的基本同态运算构建复杂算法
// 包含：幂运算、乘积、求逆、指数函数、sigmoid、同态 FFT 等
// 对应论文 §5 的扩展算法内容
// ============================================================================
class SchemeAlgo {
public:
	Scheme& scheme;  // 方案引用，提供基础同态运算接口

	SchemeAlgo(Scheme& scheme) : scheme(scheme) {};  // 构造函数：绑定方案引用


	//----------------------------------------------------------------------------------
	//   ARRAY ENCRYPTION & DECRYPTION  数组加密/解密
	//   将复数数组的每个元素分别加密为单独的密文
	//----------------------------------------------------------------------------------


	Ciphertext* encryptSingleArray(complex<double>* vals, long size);    // 复数数组 → 密文数组

	Ciphertext* encryptSingleArray(double* vals, long size);             // 实数数组 → 密文数组

	complex<double>* decryptSingleArray(SecretKey& secretKey, Ciphertext* ciphers, long size);  // 密文数组 → 复数数组


	//----------------------------------------------------------------------------------
	//   POWERS & PRODUCTS  幂运算与乘积
	//   所有幂运算都通过反复平方 + rescale 实现
	//   每次乘法后 scale 翻倍，需要 rescale 将 scale 恢复到 Δ
	//----------------------------------------------------------------------------------


	Ciphertext powerOf2(Ciphertext& cipher,  const long logDegree);      // cipher^{2^logDegree}: 反复平方 logDegree 次

	Ciphertext* powerOf2Extended(Ciphertext& cipher, const long logDegree);  // 扩展版: 返回 [cipher, cipher^2, cipher^4, ..., cipher^{2^logDegree}]

	Ciphertext power(Ciphertext& cipher, const long degree);             // cipher^degree: 递归分解为 2 的幂次 + 余项

	Ciphertext* powerExtended(Ciphertext& cipher, const long degree);    // 扩展版: 返回 [cipher^1, cipher^2, ..., cipher^degree]

	Ciphertext prodOfPo2(Ciphertext* ciphers, const long logDegree);     // ∏ ciphers (2^logDegree 个密文的乘积，树状结构)

	Ciphertext prod(Ciphertext* ciphers, const long degree);             // ∏ ciphers (任意 degree 个密文的乘积)


	//----------------------------------------------------------------------------------
	//   METHODS ON ARRAYS OF CIPHERTEXTS  密文数组操作
	//----------------------------------------------------------------------------------


	Ciphertext sum(Ciphertext* ciphers, const long size);                // 密文数组求和（逐元素同态加法）

	Ciphertext distance(Ciphertext& cipher1, Ciphertext& cipher2);       // 两个密文的距离: sum((a-b)^2) (slot 维度)

	Ciphertext* multVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size);    // 逐元素密文乘法

	void multAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size);   // 逐元素原地密文乘法

	Ciphertext* multAndModSwitchVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size);   // 逐元素乘法 + rescale

	void multModSwitchAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size);     // 逐元素原地乘法 + rescale

	Ciphertext innerProd(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size);  // 密文数组内积（先乘后加）

	Ciphertext partialSlotsSum(Ciphertext& cipher, const long slots);    // slot 求和: 将所有 slot 的值加起来（用旋转+加法）

	void partialSlotsSumAndEqual(Ciphertext& cipher, const long slots);  // 原地 slot 求和


	//----------------------------------------------------------------------------------
	//   FUNCTIONS  函数求值
	//   使用 Taylor 级数展开近似非线性函数
	//   原理: f(x) ≈ Σ a_i * x^i，先算各次幂再线性组合
	//----------------------------------------------------------------------------------


	Ciphertext inverse(Ciphertext& cipher, const long steps);            // 1/x: 基于 1/(1-(1-x)) = Σ(1-x)^i 的截断级数

	Ciphertext* inverseExtended(Ciphertext& cipher, const long steps);   // 扩展版求逆

	Ciphertext exponent(Ciphertext& cipher, long degree = 7);            // exp(x): Taylor 级数展开，默认 7 阶

	Ciphertext sigmoid(Ciphertext& cipher, long degree = 7);             // σ(x) = 1/(1+e^{-x}): Taylor 级数近似

	Ciphertext function(Ciphertext& cipher, string& funcName, const long degree);       // 通用函数求值（按名称查 Taylor 系数）

	Ciphertext* functionExtended(Ciphertext& cipher, string& funcName, const long degree);  // 扩展版通用函数求值

	//----------------------------------------------------------------------------------
	//   FFT & FFT INVERSE  同态 FFT
	//   在密文数组上执行 FFT 运算（每个元素是一个密文）
	//   使用 Cooley-Tukey 蝶形结构，用多项式乘法代替标量乘法
	//----------------------------------------------------------------------------------


	void bitReverse(Ciphertext* ciphers, const long size);               // 密文数组的 bit-reversal 置换

	void fft(Ciphertext* ciphers, const long size);                      // 同态 FFT（Cooley-Tukey 蝶形）

	void fftInvLazy(Ciphertext* ciphers, const long size);               // 同态逆 FFT（Gentleman-Sande 蝶形，不归一化）

	void fftInv(Ciphertext* ciphers, const long size);                   // 同态逆 FFT

};

#endif
