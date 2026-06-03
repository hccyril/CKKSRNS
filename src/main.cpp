/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "TestScheme.h"  // 引入测试方案类头文件
#include "Numb.h"        // 数论工具（gcd 等）

// ============================================================================
// main.cpp —— 程序入口
// 调用 TestScheme 的各种测试函数来验证 HEAAN 方案的正确性
// 每个测试函数的参数格式: (logN, L, logp, ...)
//   logN:  log₂(N)，环维度的对数
//   L:     最大 level（密文模数链长度）
//   logp:  log₂(Δ)，缩放因子的对数
//   其他参数因测试而异（slots, degree 等）
// ============================================================================

int main() {

	// 测试 1: 单值编码/解码 (logN=14, L=1, logp=55)
	// 验证 Encode + Decode 的正确性：编码一个复数，加密，解密，对比
	TestScheme::testEncodeSingle(14, 1, 55);

	// 测试 2: 批量编码/解码 (logN=15, L=6, logp=55, logSlots=3 → 8 个 slot)
	// 验证批量编码 + 加密 + 解密 + 解码的正确性
	TestScheme::testEncodeBatch(15, 6, 55, 3);

	// 测试 3: 基本同态运算 (logN=15, L=11, logp=55, logSlots=3)
	// 测试加法、乘法、常数乘法 + rescale
	TestScheme::testBasic(15, 11, 55, 3);

//	TestScheme::testConjugateBatch(15, 6, 55, 1);   // 共轭测试（已注释）

//	TestScheme::testRotateByPo2Batch(16, 26, 40, 1, 4, false);  // 2 的幂次旋转测试（已注释）

//	TestScheme::testRotateBatch(15, 6, 55, 3, 4, true);  // 任意旋转测试（已注释）

//	TestScheme::testimultBatch(16, 16, 55, 2);  // 乘虚数单位 i 测试（已注释）

//	TestScheme::testPowerOf2Batch(16, 15, 50, 2, 3);  // 2 的幂次方测试（已注释）

	// 测试 4: 求逆运算 (logN=14, L=5, logp=55, invSteps=4, logSlots=3)
	// 验证 inverse() 算法: 1/(1-x) 的截断级数展开
	TestScheme::testInverseBatch(14, 5, 55, 4, 3);

	// 测试 5: 指数运算 (logN=14, L=5, logp=55, degree=7, logSlots=3)
	// 验证 exponent() 算法: exp(x) 的 Taylor 级数展开
	TestScheme::testExponentBatch(14, 5, 55, 7, 3);

	// 测试 6: Sigmoid 运算 (logN=16, L=15, logp=55, degree=3, logSlots=3)
	// 验证 sigmoid() 算法: σ(x) = exp(x)/(1+exp(x)) 的 Taylor 近似
	TestScheme::testSigmoidBatch(16, 15, 55, 3, 3);

//	TestScheme::testSlotsSum(16, 15, 40, 3);  // slot 求和测试（已注释）

//	TestScheme::testMeanVariance(14, 3, 55, 13);  // 均值方差测试（已注释）

//	TestScheme::testHEML("data/uis.txt", 0, 5);  // 同态机器学习测试（已注释）

	return 0;  // 程序正常退出
}
