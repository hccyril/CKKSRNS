/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_CIPHERTEXT_H_  // 头文件保护宏（include guard），防止重复包含
#define HEAANNTT_CIPHERTEXT_H_  // 定义头文件保护宏名称

#include "Common.h"  // 引入公共头文件，包含 uint64_t、long 等基本类型定义

// ============================================================================
// 密文类 Ciphertext
// ============================================================================
// 对应论文 "A Full RNS Variant of Approximate Homomorphic Encryption"
// (SAC 2018, eprint 2018/931) 中的密文数据结构。
//
// 在 CKKS 方案中，一个密文 (ciphertext) 由环 R_Q = Z_Q[X]/(X^N+1) 上的
// 两个多项式 (ax, bx) 组成，即密文为 (a, b) ∈ R_Q^2。
// 在 RNS 表示下，每个多项式按其模数链 Q = q_0 * q_1 * ... * q_{l-1}
// 分解为 l 个 RNS 肢（limb），每个肢对应一个小素数模数 q_i。
// 因此每个多项式在内存中占用 N * l 个 uint64_t 空间。
//
// 密文的 level l 表示当前密文使用了多少层模数（即模数链的长度）。
// 每次同态乘法后的重线性化（relinearization）和模切换（modulus switching）
// 会消耗一层模数，使 l 减小。参见论文 Section 3.3。
// ============================================================================
class Ciphertext {

public:

	uint64_t* bx; ///< The second polynomial component b(X) of ciphertext, NTT form, RNS storage; length = N*l  // 密文的第二个多项式分量 b(X)，NTT（数论变换）表示，RNS 存储；长度 = N * l

	uint64_t* ax; ///< The first polynomial component a(X) of ciphertext, NTT form, RNS storage; length = N*l  // 密文的第一个多项式分量 a(X)，NTT（数论变换）表示，RNS 存储；长度 = N * l

	long N; ///< Dimension of Ring  // 环 R_Q = Z_Q[X]/(X^N+1) 的维度 N，N 必须是 2 的幂（论文 Section 2）

	long slots; ///< The length of plaintext vector  // 密文所编码的明文向量长度（即槽位数），slots <= N/2；对应 CKKS 编码打包的复数个数（论文 Section 3.1）

	long l; ///< The level of this ciphertext  // 密文当前所处的层级（level），等于当前模数链中剩余模数的个数；初始加密时 l = L（最大层数），每次模切换后递减（论文 Section 3.3, Figure 2）

	// Default constructor  // 默认构造函数：构造一个空的密文对象，所有指针置空、参数置零
	Ciphertext();

	// Constructor  // 参数化构造函数：用给定的多项式数组指针和密文参数直接构造密文（浅拷贝，不分配新内存）
	Ciphertext(uint64_t* ax, uint64_t* bx, long N, long slots, long l);

	// Copy constructor  // 拷贝构造函数：深拷贝另一个密文对象（分配新内存并逐元素复制所有多项式系数）
	Ciphertext(const Ciphertext& cipher);
	Ciphertext& operator=(const Ciphertext &o);  // 赋值运算符重载：深拷贝赋值，处理自赋值情况，返回当前对象引用以支持链式赋值
		
};

#endif  // 头文件保护宏结束
