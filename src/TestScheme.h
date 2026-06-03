/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_TESTSCHEME_H_     // 头文件保护宏
#define HEAANNTT_TESTSCHEME_H_

#include <iostream>                // cout 等

using namespace std;

// ============================================================================
// 测试方案类：提供 ~14 种测试函数，覆盖 HEAAN 方案的所有核心操作
// 每个测试函数都是 static 的，可直接调用
// 参数说明:
//   logN     - log₂(N)，环维度的对数
//   L        - 最大 level（密文模数链长度）
//   logp     - log₂(Δ)，缩放因子的对数
//   logSlots - log₂(slots)，slot 数的对数
// ============================================================================
class TestScheme {
public:
	// --- 基础编解码测试 ---
	static void testEncodeSingle(long logN, long L, long logp);                   // 单值编码/解码测试
	static void testEncodeBatch(long logN, long L, long logp, long logSlots);     // 批量编码/解码测试

	// --- 基本同态运算测试 ---
	static void testBasic(long logN, long L, long logp, long logSlots);           // 加法/乘法/常数乘 + rescale

	// --- 槽位操作测试 ---
	static void testConjugateBatch(long logN, long L, long logp, long logSlots);              // 共轭测试
	static void testimultBatch(long logN, long L, long logp, long logSlots);                  // 乘虚数单位 i 测试
	static void testRotateByPo2Batch(long logN, long L, long logp, long logRotSlots, long logSlots, bool isLeft);  // 2 的幂次旋转
	static void testRotateBatch(long logN, long L, long logp, long rotSlots, long logSlots, bool isLeft);          // 任意旋转
	static void testSlotsSum(long logN, long L, long logp, long logSlots);                    // slot 求和

	//----------------------------------------------------------------------------------
	//   POWER & PRODUCT TESTS  幂运算与乘积测试
	//----------------------------------------------------------------------------------

	static void testPowerOf2Batch(long logN, long L, long logp, long logDegree, long logSlots);   // 2 的幂次方
	static void testPowerBatch(long logN, long L, long logp, long degree, long logSlots);         // 任意幂次
	static void testProdOfPo2Batch(long logN, long L, long logp, long logDegree, long logSlots);  // 树状乘积
	static void testProdBatch(long logN, long L, long logp, long degree, long logSlots);          // 一般乘积

	//----------------------------------------------------------------------------------
	//   FUNCTION TESTS  函数求值测试
	//----------------------------------------------------------------------------------

	static void testInverseBatch(long logN, long L, long logp, long invSteps, long logSlots);   // 求逆
	static void testLogarithmBatch(long logN, long L, long logp, long degree, long logSlots);   // 对数
	static void testExponentBatch(long logN, long L, long logp, long degree, long logSlots);    // 指数
	static void testSigmoidBatch(long logN, long L, long logp, long degree, long logSlots);     // sigmoid

};

#endif /* TESTSCHEME_H_ */
