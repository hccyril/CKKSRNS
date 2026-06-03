/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_SECRETKEY_H_  // 头文件保护宏，防止重复包含
#define HEAANNTT_SECRETKEY_H_

#include "Common.h"            // 引入公共头文件（标准库类型定义）
#include "Context.h"           // 引入上下文类，提供采样函数 sampleHWT() 和 NTT 变换 NTTAndEqual()

// ============================================================================
// 秘密密钥类 SecretKey
// 论文 §4 KeyGen: 秘密密钥 s ∈ R 是小汉明重量的稀疏三元多项式
// 系数来自 {-1, 0, 1}，恰好有 h 个非零系数（汉明重量为 h）
// 在代码中 sx 以 NTT 形式存储，覆盖 D = B∪C 所有模数（L 个密文模数 + K 个特殊模数）
// 解密公式: Dec(sk, ct) = bx + ax·s (mod q_0)  (论文 §4 Dec)
// ============================================================================
class SecretKey {
public:

	uint64_t* sx;              // 秘密密钥多项式 s(X) 的 NTT 表示，RNS 存储
	                           // 内存大小: N * (L + K) 个 uint64_t
	                           // 前 L*N 个元素: s mod q_0, s mod q_1, ..., s mod q_{L-1}（密文模数）
	                           // 后 K*N 个元素: s mod p_0, ..., s mod p_{K-1}（特殊模数，用于 key switching）

	// 构造函数：根据 Context 参数生成秘密密钥
	// 内部调用 sampleHWT() 采样汉明重量为 h 的稀疏三元多项式
	// 然后调用 NTTAndEqual() 将所有 RNS 分量转换为 NTT 表示
	SecretKey(Context& context);

};


#endif
