/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Ciphertext.h"  // 引入密文类的头文件

// ============================================================================
// 默认构造函数
// 构造一个空的密文对象，所有指针初始化为 nullptr，所有参数初始化为 0。
// 此时密文不持有任何有效数据，需在后续通过赋值或参数化构造函数初始化。
// ============================================================================
Ciphertext::Ciphertext() : ax(nullptr), bx(nullptr), N(0), slots(0), l(0) {}  // 初始化列表：ax、bx 置空（不分配内存），N、slots、l 均置零

// ============================================================================
// 参数化构造函数
// 直接传入已有的多项式数组指针和密文参数来构造密文对象。
// 注意：此处为浅拷贝（shallow copy），即直接使用外部传入的指针，
//       不分配新内存，也不复制数据。调用者需保证指针在密文生命周期内有效。
// 参数：
//   ax    - 密文第一分量 a(X) 的多项式系数数组（NTT 域，RNS 存储），长度应为 N*l
//   bx    - 密文第二分量 b(X) 的多项式系数数组（NTT 域，RNS 存储），长度应为 N*l
//   N     - 环维度（论文 Section 2: N 为 2 的幂）
//   slots - 明文向量长度（槽位数）
//   l     - 密文层级（模数链长度）
// ============================================================================
Ciphertext::Ciphertext(uint64_t* ax, uint64_t* bx, long N, long slots, long l) : ax(ax), bx(bx), N(N), slots(slots), l(l){}  // 初始化列表：直接绑定外部传入的指针，拷贝标量参数值

// ============================================================================
// 拷贝构造函数
// 深拷贝另一个密文对象。为 ax 和 bx 分别分配新的内存（大小各为 N*l），
// 并逐元素复制所有 RNS 肢的多项式系数。
// 这保证了新密文与原密文在内存上完全独立。
// ============================================================================
Ciphertext::Ciphertext(const Ciphertext& cipher) : N(cipher.N), slots(cipher.slots), l(cipher.l) {  // 初始化列表：从源密文拷贝 N、slots、l 三个标量参数
	ax = new uint64_t[N * l];  // 为第一分量 a(X) 分配 N*l 个 uint64_t 的连续内存（N 个系数 × l 个 RNS 肢）
	bx = new uint64_t[N * l];  // 为第二分量 b(X) 分配 N*l 个 uint64_t 的连续内存
	for (long i = 0; i < N * l; ++i) {  // 遍历所有 N*l 个元素（逐位置复制，覆盖全部 RNS 肢）
		ax[i] = cipher.ax[i];  // 深拷贝第一分量 a(X) 的第 i 个元素
		bx[i] = cipher.bx[i];  // 深拷贝第二分量 b(X) 的第 i 个元素
	}
}

// ============================================================================
// 赋值运算符重载
// 实现深拷贝赋值语义。首先处理自赋值情况，然后释放当前对象已有的内存，
// 再从源对象分配新内存并逐元素复制数据。
// ============================================================================
Ciphertext& Ciphertext::operator=(const Ciphertext& o) {
	if(this == &o) return *this; // handling of self assignment, thanks for your advice, arul.  // 自赋值检测：若目标与源为同一对象则直接返回，避免后续 delete 导致数据丢失
	delete[] ax;  // 释放当前对象原有的第一分量 a(X) 数组内存，防止内存泄漏
	delete[] bx;  // 释放当前对象原有的第二分量 b(X) 数组内存
	N = o.N;  // 拷贝源对象的环维度 N
	l = o.l;  // 拷贝源对象的密文层级 l
	slots = o.slots;  // 拷贝源对象的槽位数 slots
	ax = new uint64_t[N * l];  // 为第一分量 a(X) 重新分配 N*l 大小的内存
	bx = new uint64_t[N * l];  // 为第二分量 b(X) 重新分配 N*l 大小的内存
	for (long i = 0; i < N * l; ++i) {  // 遍历所有 N*l 个元素
		ax[i] = o.ax[i];  // 深拷贝源对象第一分量 a(X) 的第 i 个元素
		bx[i] = o.bx[i];  // 深拷贝源对象第二分量 b(X) 的第 i 个元素
	}
	return *this;  // 返回当前对象的引用，支持链式赋值（如 a = b = c）
}
