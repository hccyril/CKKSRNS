/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_PLAINTEXT_H_  // 头文件保护宏（include guard），防止重复包含
#define HEAANNTT_PLAINTEXT_H_  // 定义头文件保护宏名称

#include "Common.h"  // 引入公共头文件，包含 uint64_t、long 等基本类型定义

// ============================================================================
// 明文类 Plaintext
// ============================================================================
// 对应论文 "A Full RNS Variant of Approximate Homomorphic Encryption"
// (SAC 2018, eprint 2018/931) 中的明文数据结构。
//
// 在 CKKS 方案中，明文是一个编码后的多项式 m(X) ∈ R_Q = Z_Q[X]/(X^N+1)。
// 明文多项式通过对复数向量进行 CKKS 编码得到（论文 Section 3.1, Figure 1）：
//   1. 将 slots 个复数值通过典型嵌入（canonical embedding）的逆映射编码为多项式
//   2. 乘以一个缩放因子 Delta = 2^p（论文 Section 3.2）后取整，得到整数系数多项式
//
// 在 RNS 表示下，明文多项式按模数链 Q = q_0 * q_1 * ... * q_{l-1}
// 分解为 l 个 RNS 肢（limb），每个系数对每个小素数模数 q_i 取模存储。
// 因此多项式在内存中占用 N * l 个 uint64_t 空间。
// 明文以 NTT（数论变换）形式存储，以便后续与密文进行高效的多项式运算。
// ============================================================================
class Plaintext {

public:

	uint64_t* mx;  ///< 明文多项式 m(X)，NTT（数论变换）表示，RNS 存储；长度 = N * l  // 明文多项式系数数组，NTT 域，RNS 按肢存储

	long N;  ///< 环 R_Q = Z_Q[X]/(X^N+1) 的维度 N，N 必须是 2 的幂（论文 Section 2）  // 环维度

	long slots;  ///< 明文向量的长度（槽位数），slots <= N/2；对应 CKKS 编码打包的复数个数（论文 Section 3.1）  // 明文槽位数

	long l;  ///< 明文多项式所占的模数链层数，等于 RNS 肢的个数（论文 Section 3）  // 模数层级（RNS 肢数量）

	Plaintext();  // 默认构造函数：构造空的明文对象

	Plaintext(uint64_t* mx, long N, long slots, long l);  // 参数化构造函数：用给定的多项式数组指针和参数构造明文（浅拷贝）

	Plaintext(const Plaintext& ptxt);  // 拷贝构造函数：深拷贝另一个明文对象

	Plaintext& operator=(const Plaintext &o);  // 赋值运算符重载：深拷贝赋值
		
};


#endif  // 头文件保护宏结束
