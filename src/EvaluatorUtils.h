/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAAN_EVALUATORUTILS_H_  // 头文件保护宏
#define HEAAN_EVALUATORUTILS_H_

#include "Context.h"  // 引入上下文类

// 辅助工具类：提供随机数生成和数组旋转功能
// 主要用于测试框架中生成随机测试数据、以及对比验证旋转操作的正确性
class EvaluatorUtils {
public:

	//----------------------------------------------------------------------------------
	//   RANDOM REAL AND COMPLEX NUMBERS  随机实数与复数生成
	//----------------------------------------------------------------------------------

	/**
	 * generate random real value in range (0, bound).
	 * 生成 (0, bound) 范围内的随机实数
	 */
	static double randomReal(double bound = 1.0);

	/**
	 * generate random complex value with both real and imaginary part in range (0, bound).
	 * 生成实部和虚部都在 (0, bound) 范围内的随机复数
	 */
	static complex<double> randomComplex(double bound = 1.0);

	/**
	 * generate random complex value with norm 1 and angle bound in (0, anglebound) in radiant
	 * 生成模长为 1、辐角在 (0, anglebound) 范围内的随机复数（单位圆上的点）
	 */
	static complex<double> randomCircle(double anglebound = 1.0);

	/**
	 * generate array of random real values in range (0, bound)
	 * 生成随机实数数组
	 */
	static double* randomRealArray(long size, double bound = 1.0);

	/**
	 * generate array of random complex values with both real and imaginary part in range (0, bound).
	 * 生成随机复数数组
	 */
	static complex<double>* randomComplexArray(long size, double bound = 1.0);

	/**
	 * generate array of random complex values with norm 1 and angle bound in (0, anglebound) in radiant
	 * 生成单位圆上随机复数数组
	 */
	static complex<double>* randomCircleArray(long size, double bound = 1.0);

	//----------------------------------------------------------------------------------
	//   ROTATIONS  数组旋转（用于测试验证 slot 旋转的正确性）
	//----------------------------------------------------------------------------------

	/**
	 * left indexes rotation of values
	 * @param[in, out] vals: array of values  待旋转的数组
	 * @param[in] size: array size            数组大小
	 * @param[in] rotSize: rotation size      左移位数
	 */
	static void leftRotateAndEqual(complex<double>* vals, const long size, const long rotSize);

	/**
	 * right indexes rotation of values
	 * @param[in, out] vals: array of values  待旋转的数组
	 * @param[in] size: array size            数组大小
	 * @param[in] rotSize: rotation size      右移位数
	 */
	static void rightRotateAndEqual(complex<double>* vals, const long size, const long rotSize);

};

#endif
