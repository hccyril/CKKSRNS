/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Plaintext.h"  // 引入明文类的头文件

// ============================================================================
// 默认构造函数
// 构造一个空的明文对象，所有指针初始化为 nullptr，所有参数初始化为 0。
// 此时明文不持有任何有效数据，需在后续通过赋值或参数化构造函数初始化。
// ============================================================================
Plaintext::Plaintext() : mx(nullptr), N(0), slots(0), l(0) {}  // 初始化列表：mx 置空（不分配内存），N、slots、l 均置零

// ============================================================================
// 参数化构造函数
// 直接传入已有的多项式数组指针和明文参数来构造明文对象。
// 注意：此处为浅拷贝（shallow copy），即直接使用外部传入的指针，
//       不分配新内存，也不复制数据。调用者需保证指针在明文生命周期内有效。
// 参数：
//   mx    - 明文多项式 m(X) 的系数数组（NTT 域，RNS 存储），长度应为 N*l
//   N     - 环维度（论文 Section 2: N 为 2 的幂）
//   slots - 明文向量长度（槽位数），对应 CKKS 编码的复数个数
//   l     - 模数链层数（RNS 肢的个数）
// ============================================================================
Plaintext::Plaintext(uint64_t* mx, long N, long slots, long l) : mx(mx), N(N), slots(slots), l(l) {}  // 初始化列表：直接绑定外部传入的指针，拷贝标量参数值

// ============================================================================
// 拷贝构造函数
// 深拷贝另一个明文对象。为 mx 分配新的内存（大小为 N*l），
// 并逐元素复制所有 RNS 肢的多项式系数。
// 这保证了新明文与原明文在内存上完全独立。
// ============================================================================
Plaintext::Plaintext(const Plaintext& ptxt) : N(ptxt.N), slots(ptxt.slots), l(ptxt.l){  // 初始化列表：从源明文拷贝 N、slots、l 三个标量参数
	mx = new uint64_t[N * l];  // 为明文多项式 m(X) 分配 N*l 个 uint64_t 的连续内存（N 个系数 × l 个 RNS 肢）
	for (long i = 0; i < N * l; ++i) {  // 遍历所有 N*l 个元素（逐位置复制，覆盖全部 RNS 肢）
		mx[i] = ptxt.mx[i];  // 深拷贝明文多项式 m(X) 的第 i 个元素
	}
}

// ============================================================================
// 赋值运算符重载
// 实现深拷贝赋值语义。首先处理自赋值情况，然后释放当前对象已有的内存，
// 再从源对象分配新内存并逐元素复制数据。
// ============================================================================
Plaintext& Plaintext::operator=(const Plaintext& o) {
	if(this == &o) return *this; // handling of self assignment, thanks for your advice, arul.  // 自赋值检测：若目标与源为同一对象则直接返回，避免后续 delete 导致数据丢失
	delete[] mx;  // 释放当前对象原有的明文多项式 m(X) 数组内存，防止内存泄漏
	N = o.N;  // 拷贝源对象的环维度 N
	l = o.l;  // 拷贝源对象的模数层级 l
	slots = o.slots;  // 拷贝源对象的槽位数 slots
	mx = new uint64_t[N * l];  // 为明文多项式 m(X) 重新分配 N*l 大小的内存
	for (long i = 0; i < N * l; ++i) {  // 遍历所有 N*l 个元素
		mx[i] = o.mx[i];  // 深拷贝源对象明文多项式 m(X) 的第 i 个元素
	}
	return *this;  // 返回当前对象的引用，支持链式赋值（如 a = b = c）
}
