/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_KEY_H_  // 头文件保护宏（include guard），防止重复包含
#define HEAANNTT_KEY_H_  // 定义头文件保护宏名称

#include <cstdint>  // 引入 C++ 标准整数类型头文件，提供 uint64_t 等固定宽度整数类型

using namespace std;  // 使用标准命名空间 std（注：在头文件中使用 using namespace 通常不推荐，但此处保持原有代码不变）

// ============================================================================
// 评估密钥类 Key（基类）
// ============================================================================
// 对应论文 "A Full RNS Variant of Approximate Homomorphic Encryption"
// (SAC 2018, eprint 2018/931) 中的评估密钥（evaluation key）数据结构。
//
// 在 CKKS 方案中，评估密钥用于同态乘法后的重线性化（relinearization）操作，
// 将密文从 3 个多项式分量缩减回 2 个分量。参见论文 Section 3.3。
//
// 评估密钥本质上是对密钥多项式 s 的某种编码（如 s^2 的加密），
// 由两个环多项式 (ax, bx) 组成，类似于密文结构。
// 在 RNS 表示下，每个多项式按模数链分解为多个 RNS 肢存储。
//
// 该类声明为带有虚析构函数的基类，
// 可被派生为不同类型的密钥（如重线性化密钥 rlk、旋转密钥 rot 等）。
// ============================================================================
class Key {
public:

	uint64_t* ax;  ///< 评估密钥的第一个多项式分量 a(X)，NTT 表示，RNS 存储  // 评估密钥第一分量 a(X) 的多项式系数数组

	uint64_t* bx;  ///< 评估密钥的第二个多项式分量 b(X)，NTT 表示，RNS 存储  // 评估密钥第二分量 b(X) 的多项式系数数组

	Key(uint64_t* ax, uint64_t* bx);  // 构造函数：用给定的两个多项式数组指针构造评估密钥（浅拷贝）

	virtual ~Key();  // 虚析构函数：声明为 virtual 以支持派生类的正确多态析构（确保 delete 基类指针时能调用派生类析构函数）
};

#endif /* KEY_H_ */  // 头文件保护宏结束
