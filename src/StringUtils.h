/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAAN_STRINGUTILS_H_  // 头文件保护宏
#define HEAAN_STRINGUTILS_H_

#include <complex>   // complex<double> 复数类型

#include "Common.h"  // 公共头文件

using namespace std;

// ============================================================================
// 字符串工具类：提供数组打印和结果对比功能
// 所有方法均为 static，用于测试框架中打印明文/密文对比结果
// 格式: m前缀 = 明文(message), d前缀 = 解密结果(decrypted), e前缀 = 误差(error)
// ============================================================================
class StringUtils {
public:


	//----------------------------------------------------------------------------------
	//   SHOW ARRAY  数组打印
	//----------------------------------------------------------------------------------


	static void show(uint64_t* vals, long size);  // 打印 uint64_t 数组
	/**
	 * prints in console array
	 * @param[in] vals: long array
	 * @param[in] size: array size
	 */
	static void show(long* vals, long size);      // 打印 long 数组

	/**
	 * prints in console array
	 * @param[in] vals: double array
	 * @param[in] size: array size
	 */
	static void show(double* vals, long size);    // 打印 double 数组

	/**
	 * prints in console array
	 * @param[in] vals: complex array
	 * @param[in] size: array size
	 */
	static void show(complex<double>* vals, long size);  // 打印复数数组


	//----------------------------------------------------------------------------------
	//   SHOW & COMPARE ARRAY  打印并对比数组（明文 vs 解密结果）
	//----------------------------------------------------------------------------------


	/**
	 * prints in console val1, val2 and (val1-val2)
	 * 打印两个值及其差（误差）
	 * @param[in] val1: double value   明文值
	 * @param[in] val2: double value   解密值
	 * @param[in] prefix: string prefix  前缀标识
	 */
	static void showcompare(double val1, double val2, string prefix);

	/**
	 * prints in console val1, val2 and (val1-val2)
	 * @param[in] val1: complex value  明文复数
	 * @param[in] val2: complex value  解密复数
	 * @param[in] prefix: string prefix  前缀标识
	 */
	static void showcompare(complex<double> val1, complex<double> val2, string prefix);

	/**
	 * prints in console pairwise val1[i], val2[i] and (val1[i]-val2[i])
	 * 逐元素打印两个数组及其差
	 * @param[in] vals1: double array   明文数组
	 * @param[in] vals2: double array   解密数组
	 * @param[in] size: array size      数组大小
	 * @param[in] prefix: string prefix  前缀标识
	 */
	static void showcompare(double* vals1, double* vals2, long size, string prefix);

	/**
	 * prints in console pairwise val1[i], val2[i] and (val1[i]-val2[i])
	 * 逐元素打印两个复数数组及其差
	 * @param[in] vals1: complex array  明文复数数组
	 * @param[in] vals2: complex array  解密复数数组
	 * @param[in] size: array size
	 * @param[in] prefix: string prefix
	 */
	static void showcompare(complex<double>* vals1, complex<double>* vals2, long size, string prefix);

	// 以下 4 个重载用于标量与数组的对比
	static void showcompare(double* vals1, double val2, long size, string prefix);        // 数组 vs 标量
	static void showcompare(complex<double>* vals1, complex<double> val2, long size, string prefix);  // 复数数组 vs 标量
	static void showcompare(double val1, double* vals2, long size, string prefix);        // 标量 vs 数组
	static void showcompare(complex<double> val1, complex<double>* vals2, long size, string prefix);  // 标量 vs 复数数组

};

#endif
