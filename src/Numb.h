/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_NUMB_H_       // 头文件保护宏，防止重复包含
#define HEAANNTT_NUMB_H_

// 数论工具库：提供底层模运算原语
// 这些函数是整个 HEAAN 方案的数学基础设施
// 所有多项式环 R_q 中的运算最终都调用这些函数

#include <iostream>             // cout 等
#include <vector>               // 容器
#include <set>                  // 集合（用于找素因子）
#include <math.h>               // sqrt 等数学函数
#include <stdint.h>             // uint64_t 等固定宽度整数类型

#include "Common.h"             // 公共头文件

using namespace std;

// ---------- 基本模运算 ----------

void negate(uint64_t& r, uint64_t a);                          // r = -a（取负，等价于 p - a）

void addMod(uint64_t& r, uint64_t a, uint64_t b, uint64_t p); // r = (a + b) mod p  模加法

void addModAndEqual(uint64_t& a, uint64_t b, uint64_t p);     // a = (a + b) mod p  原地模加

void subMod(uint64_t& r, uint64_t a, uint64_t b, uint64_t p); // r = (a - b) mod p  模减法

void subModAndEqual(uint64_t& a, uint64_t b, uint64_t p);     // a = (a - b) mod p  原地模减

void mulMod(uint64_t& r, uint64_t a, uint64_t b, uint64_t p); // r = (a * b) mod p  模乘法（使用 __int128 防溢出）

// Barrett reduction 加速的模乘：避免大整数除法
// pr = Barrett 预计算常数 ⌊2^{twok} / p⌋
// twok = 2 * (⌊log₂(p)⌋ + 1)，位移量
void mulModBarrett(uint64_t& r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr, long twok);

// Barrett reduction 模运算（64 位输入）
void modBarrett(uint64_t &r, uint64_t a, uint64_t m, uint64_t mr, long twok);

// Barrett reduction 模运算（128 位输入，用于 NTT 蝶形中的大中间值）
void modBarrett(uint64_t &r, unsigned __int128 a, uint64_t m, uint64_t mr, long twok);

void mulModAndEqual(uint64_t& a, uint64_t b, uint64_t p);     // a = (a * b) mod p  原地模乘

// ---------- 模逆与模幂 ----------

uint64_t invMod(uint64_t x, uint64_t p);   // x^{-1} mod p  模逆（Fermat 小定理: x^{p-2} mod p）
                                            // 论文 §2.2: 基转换中的 qhat_j^{-1} mod q_j

uint64_t powMod(uint64_t x, uint64_t y, uint64_t modulus); // x^y mod modulus  快速模幂（平方-乘法）
                                                             // 论文 §4 KeyGen 中的模幂运算

uint64_t inv(uint64_t x);     // x^{-1} mod 2^64  模逆（用于 Montgomery reduction）

uint64_t pow(uint64_t x, uint64_t y);  // x^y  普通整数快速幂（无模）

// ---------- NTT 相关 ----------

uint32_t bitReverse(uint32_t x);  // 32 位整数的位反转
                                   // 论文 §2.2 NTT: Cooley-Tukey 蝶形运算需要对输入做 bit-reversal 置换

// ---------- 数论基础 ----------

uint64_t gcd(uint64_t a, uint64_t b);  // 最大公约数（Euclidean 算法）

long gcd(long a, long b);              // 最大公约数（long 重载）

void findPrimeFactors(set<uint64_t> &s, uint64_t number);  // 求一个数的所有素因子（存入集合 s）
                                                             // 用于 findPrimitiveRoot 中的原根判定

uint64_t findPrimitiveRoot(uint64_t modulus);  // 找模 modulus 的原根（生成元）
                                                // 论文 §2.2 NTT: 需要 M 次本原单位根
                                                // 原根 g 满足: g^k mod p 互不相同 (k=1,...,p-1)

uint64_t findMthRootOfUnity(uint64_t M, uint64_t mod);  // 找模 mod 的 M 次本原单位根
                                                          // 论文 §2.2 NTT 旋转因子: ω^M ≡ 1 (mod p)
                                                          // 计算: ω = g^((p-1)/M) mod p，其中 g 是原根

bool primeTest(uint64_t p);  // Miller-Rabin 素数检测（200 轮）
                              // 论文 §4 Setup: 生成密文模数 q_j 和特殊模数 p_i 时需要验证素性
                              // 所有模数必须满足 q ≡ 1 (mod 2N) 以支持 NTT

#endif
