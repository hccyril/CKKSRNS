/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "SecretKey.h"  // 引入秘密密钥类的头文件

// ============================================================================
// 秘密密钥构造函数 —— 论文 §4 KeyGen 的第一步: s ← χ_key
// ============================================================================
// 生成流程（对应论文 §4 Setup + KeyGen）:
//   1. 分配内存: sx 需要在全部 D = B∪C 模数上存储（L 个密文模数 + K 个特殊模数）
//   2. 采样: 用汉明重量采样 (HWT) 生成稀疏三元多项式
//      - 系数来自 {-1, 0, 1}
//      - 恰好有 h 个非零系数（默认 h = 64，见 Context 构造函数参数）
//      - 论文附录 A 噪声分析: 小汉明重量使 ||s|| 可控
//   3. NTT: 将所有 RNS 分量（L+K 个模数）从系数表示转为 NTT 表示
//      后续所有涉及 s 的运算（解密、key switching）都在 NTT 域进行
// ============================================================================
SecretKey::SecretKey(Context& context) {
	// 分配 sx 数组: N * (L + K) 个 uint64_t，初始化为 0
	// L = 密文模数个数, K = 特殊模数个数
	// () 语法保证零初始化，未赋值的元素为 0
	sx = new uint64_t[context.N * (context.L + context.K)]();  // 论文 §4: sx 在 D = B∪C 基上存储

	// 论文 §4 χ_key 分布: 汉明重量为 h 的稀疏三元多项式采样
	// 参数: 输出指针 sx, 密文模数层数 L, 特殊模数个数 K
	// 内部实现: 随机选 h 个位置，每个位置以等概率赋 +1 或 -1
	context.sampleHWT(sx, context.L, context.K);  // 论文 §4: s ← χ_key (汉明重量采样)

	// 将 sx 的所有 RNS 分量（L 个密文模数 + K 个特殊模数）转为 NTT 表示
	// NTT 表示下的多项式乘法是逐点相乘，复杂度 O(N)
	// 若保持系数表示则乘法需要 O(N²) 或 O(N log N) 的 FFT
	context.NTTAndEqual(sx, context.L, context.K);  // 论文 §2.2: 所有分量转 NTT 表示
}
