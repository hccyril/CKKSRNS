/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "Numb.h"  // 引入数论工具库头文件

// ============================================================================
// 取负运算
// ============================================================================
void negate(uint64_t &r, uint64_t a) {
	r = -a;  // 对 uint64_t 取负等价于 2^64 - a，即模 2^64 下的取负
}

// ============================================================================
// 模加法: r = (a + b) mod m
// ============================================================================
void addMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t m) {
	r = (a + b) % m;  // 直接相加后取模（当 a,b < m 时，a+b < 2m，不会溢出 uint64_t）
}

// ============================================================================
// 模减法: r = (a - b) mod m
// 技巧: 先对 b 取模，再用 a + (m - b%m) 避免负数
// ============================================================================
void subMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t m) {
	r = b % m;              // 先对 b 取模，保证 b < m
	r = (a + m - r) % m;    // a + m - b 保证非负，再取模
	// 等价于 (a - b) mod m，但避免了 uint64_t 的负数问题
}

// ============================================================================
// 模乘法: r = (a * b) mod m
// 使用 unsigned __int128 避免中间乘积溢出 uint64_t
// a * b 最大可达 (2^64-1)^2 ≈ 2^128，需要 128 位整数
// ============================================================================
void mulMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t m) {
	unsigned __int128 mul = static_cast<unsigned __int128>(a % m) * (b % m);  // 先取模再乘，用 128 位存储中间结果
	mul %= static_cast<unsigned __int128>(m);                                  // 对 m 取模
	r = static_cast<uint64_t>(mul);                                            // 截断回 64 位
}

// ============================================================================
// Barrett Reduction 加速的模乘法
// Barrett reduction 原理: 用预计算的 μ = ⌊2^{twok}/p⌋ 替代除法
//   商 q ≈ ⌊a * μ / 2^{twok}⌋ （用移位代替除法）
//   余数 r = a - q * p
//   若 r ≥ p 则修正: r -= p
// 参数:
//   a, b: 输入操作数
//   p:    模数
//   pr:   Barrett 预计算常数 μ = ⌊2^{twok}/p⌋
//   twok: 位移量 = 2*(⌊log₂(p)⌋+1)
// ============================================================================
void mulModBarrett(uint64_t& r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr, long twok) {
	unsigned __int128 mul = static_cast<unsigned __int128>(a % p) * (b % p);  // 128 位中间乘积
	modBarrett(r, mul, p, pr, twok);                                           // 调用 Barrett reduction
}

// ============================================================================
// Barrett Reduction（64 位输入版本）
// ============================================================================
void modBarrett(uint64_t &r, uint64_t a, uint64_t m, uint64_t mr, long twok) {
	unsigned __int128 tmp = static_cast<unsigned __int128>(a) * mr;  // 计算 a * μ（128 位）
	tmp >>= twok;                                                     // 右移 twok 位，得到商的近似值 q
	tmp *= m;                                                         // q * m
	tmp = a - tmp;                                                    // 余数 r = a - q*m
	r = static_cast<uint64_t>(tmp);                                   // 截断为 64 位
	if(r < m) {                                                       // 如果 r < m，说明商正确
		return;
	}
	else {                                                            // 否则商偏小 1，修正
		r -= m;
		return;
	}

}

// ============================================================================
// Barrett Reduction（128 位输入版本）
// 用于 NTT 蝶形运算中，中间值是 128 位的
// 计算逻辑与 64 位版本相同，但需要处理 128 位的拆分
// ============================================================================
void modBarrett(uint64_t &r, unsigned __int128 a, uint64_t m, uint64_t mr, long twok) {
	uint64_t atop, abot;
	abot = static_cast<uint64_t>(a);              // 取 a 的低 64 位
	atop = static_cast<uint64_t>(a >> 64);        // 取 a 的高 64 位
	// 计算 a * mr:
	// 先算 abot * mr，取高 64 位（相当于右移 64 位）
	unsigned __int128 tmp = static_cast<unsigned __int128>(abot) * mr;
	tmp >>= 64;
	// 再加上 atop * mr
	tmp += static_cast<unsigned __int128>(atop) * mr;
	tmp >>= twok - 64;  // 再右移 twok-64 位，总共右移 twok 位
	tmp *= m;            // q * m
	tmp = a - tmp;       // 余数 r = a - q*m
	r = static_cast<uint64_t>(tmp);
	if(r >= m) r-= m;    // 修正：如果 r ≥ m，则 r -= m
}

// ============================================================================
// 模逆: x^{-1} mod m
// 使用 Fermat 小定理: 当 m 是素数时，x^{-1} = x^{m-2} mod m
// 前提: gcd(x, m) = 1，否则逆元不存在
// 论文 §2.2: 基转换中需要 qhat_j^{-1} mod q_j
// ============================================================================
uint64_t invMod(uint64_t x, uint64_t m) {
	uint64_t temp = x % m;
	if (gcd(temp, m) != 1) {                                       // 检查 x 与 m 是否互素
		throw invalid_argument("Inverse doesn't exist!!!");        // 不互素则无逆元，抛异常
	} else {
		return powMod(temp, m - 2, m);                             // x^{m-2} mod m = x^{-1} mod m
	}
}

// ============================================================================
// 快速模幂: x^y mod modulus
// 使用平方-乘法（square-and-multiply）算法，时间复杂度 O(log y)
// 算法: 将 y 表示为二进制，从高位到低位逐位处理
//   每一位: res = res^2 * x^{y_i} mod modulus
// ============================================================================
uint64_t powMod(uint64_t x, uint64_t y, uint64_t modulus) {
	uint64_t res = 1;                    // 初始化结果
	while (y > 0) {
		if (y & 1) {                     // 若 y 的最低位为 1
			mulMod(res, res, x, modulus); // res = res * x mod modulus
		}
		y = y >> 1;                      // y 右移一位（除以 2）
		mulMod(x, x, x, modulus);        // x = x^2 mod modulus
	}
	return res;
}

// ============================================================================
// 模 2^64 的逆元（无模参数版本）
// 用于 Montgomery reduction
// ============================================================================
uint64_t inv(uint64_t x) {
	return pow(x, static_cast<uint64_t>(-1));  // -1 转为 uint64_t 是 2^64-1，即求 x^{2^64-1}
}

// ============================================================================
// 普通整数快速幂: x^y（无模）
// 同样使用平方-乘法算法
// ============================================================================
uint64_t pow(uint64_t x, uint64_t y) {
	uint64_t res = 1;
	while (y > 0) {
		if (y & 1) {        // y 的最低位为 1
			res *= x;        // res = res * x
		}
		y = y >> 1;          // y 右移一位
		x *= x;              // x = x^2
	}
	return res;
}

// ============================================================================
// 32 位整数的位反转 (bit-reversal)
// 使用分治法逐层交换位:
//   第 1 层: 相邻 1 位交换 (0xaaaaaaaa / 0x55555555)
//   第 2 层: 相邻 2 位交换 (0xcccccccc / 0x33333333)
//   第 3 层: 相邻 4 位交换 (0xf0f0f0f0 / 0x0f0f0f0f)
//   第 4 层: 相邻 8 位交换 (0xff00ff00 / 0x00ff00ff)
//   第 5 层: 高 16 位与低 16 位交换
// 论文 §2.2 NTT: Cooley-Tukey 蝶形运算需要将输入按 bit-reversal 重排
// ============================================================================
uint32_t bitReverse(uint32_t x) {
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));   // 交换相邻 1 位
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));   // 交换相邻 2 位
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));   // 交换相邻 4 位
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));   // 交换相邻 8 位
	return ((x >> 16) | (x << 16));                              // 交换高 16 位和低 16 位
}

// ============================================================================
// 最大公约数 (GCD) - Euclidean 算法
// gcd(a, b) = gcd(b%a, a)，递归直到 a = 0
// ============================================================================
uint64_t gcd(uint64_t a, uint64_t b) {
	if (a == 0) {
		return b;            // a = 0 时，gcd(0, b) = b
	}
	return gcd(b % a, a);    // 递归: gcd(a, b) = gcd(b%a, a)
}

// long 版本的重载
long gcd(long a, long b) {
	if (a == 0) {
		return b;
	}
	return gcd(b % a, a);
}

// ============================================================================
// 求一个数的所有素因子
// 将 number 的所有素因子存入集合 s
// 用于 findPrimitiveRoot 中的原根判定
// ============================================================================
void findPrimeFactors(set<uint64_t> &s, uint64_t number) {
	while (number % 2 == 0) {    // 先提取所有因子 2
		s.insert(2);
		number /= 2;
	}
	for (uint64_t i = 3; i < sqrt(number); i++) {  // 从 3 开始，步长为 2（奇数）
		while (number % i == 0) {
			s.insert(i);
			number /= i;
		}
	}
	if (number > 2) {            // 若剩余的 number > 2，则它本身是一个素因子
		s.insert(number);
	}
}

// ============================================================================
// 找模 modulus 的原根（primitive root / 生成元）
// 原根 g 满足: {g^1, g^2, ..., g^{φ(m)}} = {1, 2, ..., m-1} (mod m)
// 等价条件: 对 φ(m) 的每个素因子 q，g^{φ(m)/q} ≠ 1 (mod m)
// 论文 §2.2 NTT: 需要 M 次本原单位根，先找原根再求幂
// ============================================================================
uint64_t findPrimitiveRoot(uint64_t modulus) {
	set<uint64_t> s;
	uint64_t phi = modulus - 1;               // 对素数 modulus，φ(modulus) = modulus - 1
	findPrimeFactors(s, phi);                 // 求 φ(modulus) 的所有素因子
	for (uint64_t r = 2; r <= phi; r++) {     // 从 2 开始逐个尝试
		bool flag = false;
		for (auto it = s.begin(); it != s.end(); it++) {
			if (powMod(r, phi / (*it), modulus) == 1) {  // 若 r^{φ/q} ≡ 1 (mod m)，则 r 不是原根
				flag = true;
				break;
			}
		}
		if (flag == false) {                  // 若对所有素因子 q 都 ≠ 1，则 r 是原根
			return r;
		}
	}
	return -1;                                // 理论上不会到这里（素数一定有原根）
}

// ============================================================================
// 找模 mod 的 M 次本原单位根
// ω 满足: ω^M ≡ 1 (mod p) 且 ω^k ≢ 1 (k < M)
// 计算: 先找原根 g，然后 ω = g^((p-1)/M) mod p
// 前提: (p-1) % M == 0，否则不存在 M 次单位根
// 论文 §2.2 NTT: 蝶形运算需要 M 次本原单位根作为旋转因子
// ============================================================================
uint64_t findMthRootOfUnity(uint64_t M, uint64_t mod) {
    uint64_t res;
    res = findPrimitiveRoot(mod);                  // 先找模 mod 的原根 g
    if((mod - 1) % M == 0) {                       // 检查 (mod-1) 是否能被 M 整除
        uint64_t factor = (mod - 1) / M;           // 计算 (mod-1)/M
        res = powMod(res, factor, mod);            // ω = g^{(mod-1)/M} mod mod
        return res;
    }
    else {
        return -1;                                 // 不存在 M 次单位根
    }
}

// ============================================================================
// Miller-Rabin 素数检测
// 概率性素数检测算法，200 轮检测
// 原理: 对素数 p，将 p-1 分解为 2^s * d（d 为奇数）
//   对随机 a ∈ [1, p-1]:
//     若 a^d ≡ 1 (mod p) 或 a^{2^r * d} ≡ -1 (mod p) (某 r)，则 p "可能是素数"
//   200 轮检测均通过则 p 是素数的概率极高（误判率 < 4^{-200}）
// 论文 §4 Setup: 生成密文模数 q_j 和特殊模数 p_i 时使用
// ============================================================================
bool primeTest(uint64_t p) {
	if(p < 2) return false;                        // 小于 2 的数不是素数
	if(p != 2 && p % 2 == 0) return false;         // 除 2 外的偶数不是素数
	uint64_t s = p - 1;
	while(s % 2 == 0) {                            // 将 p-1 分解为 2^s * d
		s /= 2;
	}
	for(long i = 0; i < 200; i++) {                // 200 轮检测
		uint64_t temp1 = rand();
		temp1  = (temp1 << 32) | rand();           // 生成 64 位随机数
		temp1 = temp1 % (p - 1) + 1;              // temp1 ∈ [1, p-1]
		uint64_t temp2 = s;
		uint64_t mod = powMod(temp1,temp2,p);      // mod = temp1^s mod p
		while (temp2 != p - 1 && mod != 1 && mod != p - 1) {
			mulMod(mod, mod, mod, p);              // mod = mod^2 mod p
		    temp2 *= 2;                             // temp2 *= 2
		}
		if (mod != p - 1 && temp2 % 2 == 0) return false;  // 检测失败，p 是合数
	}
	return true;                                   // 200 轮均通过，p 极大概率是素数
}
