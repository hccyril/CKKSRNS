/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_SCHEME_H_   // 头文件保护宏
#define HEAANNTT_SCHEME_H_

#include <map>               // map 容器（用于存储密钥对）
#include <chrono>            // 时间相关

#include "Common.h"          // 公共头文件
#include "Ciphertext.h"      // 密文类
#include "Context.h"         // 全局上下文类
#include "Plaintext.h"       // 明文类
#include "SecretKey.h"       // 秘密密钥类
#include "Key.h"             // 评估密钥类
#include "Numb.h"            // 数论工具

using namespace std;

// 密钥类型标识符（作为 keyMap 的键）
static long ENCRYPTION = 0;       // 公钥加密密钥 pk
static long MULTIPLICATION  = 1;  // 乘法重线性化密钥 evk_mult (KSGen(s²,s))
static long CONJUGATION = 2;      // 共轭密钥 evk_conj (KSGen(s̄,s))

// ============================================================================
// 方案核心类 Scheme —— 论文 §4 的完整实现
// ============================================================================
// 承载 HEAAN/CKKS 方案的所有核心算法:
//   - 密钥生成: addEncKey (§4 pk), addMultKey (§4 KSGen), addConjKey, addLeftRotKey
//   - 加密解密: encryptMsg, decryptMsg (§4 Enc/Dec)
//   - 同态运算: add, sub, mult, square (§4 Add/Sub/Mult)
//   - 常数操作: addConst, multByConst (§4 常数运算)
//   - 重缩放:   reScaleBy, reScaleTo (§4 RS)
//   - 模切换:   modDownBy, modDownTo (§3.2 ModDown)
//   - 槽位操作: leftRotate, conjugate (§4 旋转/共轭)
// ============================================================================
class Scheme {
public:

	Context& context;  // 全局上下文引用，提供环运算、NTT、采样等基础功能

	map<long, Key> keyMap; ///< contain Encryption, Multiplication and Conjugation keys, if generated  密钥映射表：存储公钥/乘法密钥/共轭密钥
	map<long, Key> leftRotKeyMap; ///< contain left rotation keys, if generated  左旋转密钥映射表：rot → evk_rot

	// ========== 构造函数 ==========
	Scheme(Context& context);  // 仅绑定上下文（不自动生成密钥）

	Scheme(SecretKey& secretKey, Context& context);  // 绑定上下文 + 自动生成加密密钥和乘法密钥

	/**
	 * generates key for public encryption (key is stored in keyMap)
	 * 生成公钥加密密钥 pk（论文 §4 KeyGen 的 pk 部分）
	 * pk = (b, a), b = -a·s + e
	 */
	void addEncKey(SecretKey& secretKey);

	/**
	 * generates key for multiplication (key is stored in keyMap)
	 * 生成乘法重线性化密钥 evk（论文 §4 KSGen(s²,s)）
	 * evk = (b, a), b + a·s ≈ P·s² + e
	 */
	void addMultKey(SecretKey& secretKey);

	/**
	 * generates key for conjugation (key is stored in keyMap)
	 * 生成共轭密钥（论文 §4 KSGen(s̄,s)）
	 * evk_conj = (b, a), b + a·s ≈ P·conjugate(s) + e
	 */
	void addConjKey(SecretKey& secretKey);

	/**
	 * generates key for left rotation (key is stored in leftRotKeyMap)
	 * 生成左旋转密钥（论文 §4 KSGen(s_rot,s)）
	 * evk_rot = (b, a), b + a·s ≈ P·rotate(s,rot) + e
	 */
	void addLeftRotKey(SecretKey& secretKey, long rot);

	/**
	 * generates all keys for power-of-two left rotations (keys are stored in leftRotKeyMap)
	 * 生成所有 2 的幂次左旋转密钥: rot = 1, 2, 4, ..., N/4
	 */
	void addLeftRotKeys(SecretKey& secretKey);

	/**
	 * generates all keys for power-of-two right rotations (keys are stored in leftRotKeyMap)
	 * 生成所有 2 的幂次右旋转密钥: rot = N/2-1, N/2-2, N/2-4, ...
	 */
	void addRightRotKeys(SecretKey& secretKey);

	// ========== 编码 / 解码（论文 §2.1）==========
	Plaintext encode(double* vals, long slots, long l);             // 实数数组编码
	Plaintext encode(complex<double>* vals, long slots, long l);    // 复数数组编码
	Plaintext encodeSingle(complex<double> val, long l);            // 单复数编码

	complex<double>* decode(Plaintext& msg);                        // 批量解码
	complex<double> decodeSingle(Plaintext& msg);                   // 单值解码

	// ========== 加密 / 解密（论文 §4 Enc / Dec）==========
	// Encryption (secret and public version)
	Ciphertext encryptMsg(SecretKey& secretkey, Plaintext& message);  // 秘密密钥加密（未实现，返回空）

	Ciphertext encryptMsg(Plaintext& message);                        // 公钥加密（论文 §4 Enc_pk）

	Plaintext decryptMsg(SecretKey& secretkey, Ciphertext& cipher);   // 解密（论文 §4 Dec_sk）

	Ciphertext encrypt(double* vals, long slots, long l);             // 实数加密（编码+加密）
	Ciphertext encrypt(complex<double>* vals, long slots, long l);    // 复数加密
	Ciphertext encryptSingle(complex<double> val, long l);            // 单值加密

	complex<double>* decrypt(SecretKey& secretKey, Ciphertext& cipher);  // 批量解密（解密+解码）
	complex<double> decryptSingle(SecretKey& secretKey, Ciphertext& cipher);  // 单值解密

	// ========== 同态取反（论文 §4 Negate）==========
	// Homomorphic Negation
	Ciphertext negate(Ciphertext& cipher);         // ct_neg = (-c0, -c1)
	void negateAndEqual(Ciphertext& cipher);       // 原地取反

	// ========== 同态加法（论文 §4 Add）==========
	// Homomorphic Addition
	Ciphertext add(Ciphertext& cipher1, Ciphertext& cipher2);         // ct_add = ct1 + ct2
	void addAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);       // ct1 += ct2

	// ========== 同态减法（论文 §4 Sub）==========
	// Homomorphic Substraction
	Ciphertext sub(Ciphertext& cipher1, Ciphertext& cipher2);         // ct_sub = ct1 - ct2

	void subAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);       // ct1 -= ct2
	void sub2AndEqual(Ciphertext& cipher1, Ciphertext& cipher2);      // ct2 = ct1 - ct2

	// ========== 同态乘法（论文 §4 Mult）==========
	// Homomorphic Multiplication
	Ciphertext mult(Ciphertext& cipher1, Ciphertext& cipher2);        // ct_mult = ct1 ⊗ ct2 (含重线性化)
	void multAndEqual(Ciphertext& cipher1, Ciphertext& cipher2);      // ct1 ⊗= ct2

	// ========== 同态平方（论文 §4 Mult 的特例）==========
	// Homomorphic Squaring
	Ciphertext square(Ciphertext& cipher);         // ct_sq = ct ⊗ ct (优化版平方)
	void squareAndEqual(Ciphertext& cipher);       // 原地平方

	// ========== 常数同态运算（论文 §4 常数操作）==========
	// Homomorphic Operations with Constant

	Ciphertext imult(Ciphertext& cipher);          // 乘虚数单位 i: a(X) → a(X)·X^{N/2}
	void imultAndEqual(Ciphertext& cipher);

	Ciphertext idiv(Ciphertext& cipher);           // 除以 i（未实现）
	void idivAndEqual(Ciphertext& cipher);

	Ciphertext addConst(Ciphertext& cipher, double cnst);                   // 加实数常数
	Ciphertext addConst(Ciphertext& cipher, complex<double> cnst);          // 加复数常数（未实现）

	void addConstAndEqual(Ciphertext& cipher, double cnst);                 // 原地加实数常数
	void addConstAndEqual(Ciphertext& cipher, complex<double> cnst);        // 原地加复数常数（未实现）

	void addPcAndEqual(Ciphertext& cipher);   // 加预计算常数 pccoeff
	void addP2AndEqual(Ciphertext& cipher);   // 加 Δ²（用于 exponent 中的常数项 1）
	void addP2hAndEqual(Ciphertext& cipher);  // 加 Δ²/2（用于 sigmoid 中的常数项 1/2）

	Ciphertext multByConst(Ciphertext& cipher, double cnst);                   // 乘实数常数
	Ciphertext multByConst(Ciphertext& cipher, complex<double> cnst);          // 乘复数常数

	Ciphertext multByConstVec(Ciphertext& cipher, double* cnstVec, long slots);      // 乘常数向量（编码后逐点乘）
	Ciphertext multByConstVec(Ciphertext& cipher, complex<double>* cnstVec, long slots);

	void multByConstAndEqual(Ciphertext& cipher, double cnst);                       // 原地乘常数
	void multByConstAndEqual(Ciphertext& cipher, complex<double> cnst);

	void multByPolyAndEqual(Ciphertext& cipher, uint64_t* poly);  // 乘明文多项式

	Ciphertext multByMonomial(Ciphertext& cipher, long mdeg);     // 乘单项式 X^{mdeg}
	void multByMonomialAndEqual(Ciphertext& cipher, long mdeg);

	// ========== 重缩放（论文 §4 RS）==========
	Ciphertext reScaleBy(Ciphertext& cipher, long dl);         // 降 dl 层 rescale
	void reScaleByAndEqual(Ciphertext& cipher, long dl);       // 原地降 dl 层

	Ciphertext reScaleTo(Ciphertext& cipher, long l);          // rescale 到 level l
	void reScaleToAndEqual(Ciphertext& cipher, long l);

	// ========== 模切换（论文 §3.2 ModDown：仅丢弃模数，不做除法）==========
	Ciphertext modDownBy(Ciphertext& cipher, long dl);         // 丢弃最高 dl 个模数
	void modDownByAndEqual(Ciphertext& cipher, long dl);

	Ciphertext modDownTo(Ciphertext& cipher, long dl);         // 切换到 level dl
	void modDownToAndEqual(Ciphertext& cipher, long dl);

	// ========== 槽位旋转（论文 §4 LeftRotate / RightRotate）==========
	Ciphertext leftRotateFast(Ciphertext& cipher, long rotSlots);         // 快速左旋转（含 key switching）
	void leftRotateAndEqualFast(Ciphertext& cipher, long rotSlots);

	Ciphertext leftRotateByPo2(Ciphertext& cipher, long logRotSlots);     // 2 的幂次左旋转
	void leftRotateByPo2AndEqual(Ciphertext& cipher, long logRotSlots);

	Ciphertext rightRotateByPo2(Ciphertext& cipher, long logRotSlots);    // 2 的幂次右旋转
	void rightRotateByPo2AndEqual(Ciphertext& cipher, long logRotSlots);

	Ciphertext leftRotate(Ciphertext& cipher, long rotSlots);             // 一般左旋转（分解为 2 的幂次）
	void leftRotateAndEqual(Ciphertext& cipher, long rotSlots);

	Ciphertext rightRotate(Ciphertext& cipher, long rotSlots);            // 一般右旋转
	void rightRotateAndEqual(Ciphertext& cipher, long rotSlots);

	// ========== 共轭（论文 §4 Conjugate）==========
	Ciphertext conjugate(Ciphertext& cipher);         // 复共轭: 每个 slot 取 conjugate
	void conjugateAndEqual(Ciphertext& cipher);

};

#endif
