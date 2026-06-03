/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include <cmath>       // cos, sin, M_PI 等数学函数
#include <complex>     // complex<double> 复数类型
#include <cstdlib>     // rand() 随机数

#include "EvaluatorUtils.h"  // 引入本类头文件

//----------------------------------------------------------------------------------
//   RANDOM REAL AND COMPLEX NUMBERS  随机数生成函数实现
//----------------------------------------------------------------------------------

// 生成 (0, bound) 范围内的随机实数
// rand()/(RAND_MAX) 产生 (0,1) 的均匀分布随机数
double EvaluatorUtils::randomReal(double bound)  {
	return (double) rand()/(RAND_MAX) * bound;  // 线性映射到 (0, bound)
}

// 生成实部和虚部都在 (0, bound) 范围内的随机复数
complex<double> EvaluatorUtils::randomComplex(double bound) {
	complex<double> res;
	res.real(randomReal(bound));  // 实部 = (0, bound) 的随机实数
	res.imag(randomReal(bound));  // 虚部 = (0, bound) 的随机实数
	return res;
}

// 生成单位圆上的随机复数: cos(θ) + i·sin(θ)，θ ∈ (0, anglebound)
complex<double> EvaluatorUtils::randomCircle(double anglebound) {
	double angle = randomReal(anglebound);        // 随机辐角 θ ∈ (0, anglebound)
	complex<double> res;
	res.real(cos(angle * 2 * M_PI));              // 实部 = cos(2π·θ)
	res.imag(sin(angle * 2 * M_PI));              // 虚部 = sin(2π·θ)
	return res;
}

// 生成随机实数数组
double* EvaluatorUtils::randomRealArray(long size, double bound) {
	double* res = new double[size];               // 分配数组内存
	for (long i = 0; i < size; ++i) {
		res[i] = randomReal(bound);               // 逐个生成随机实数
	}
	return res;
}

// 生成随机复数数组
complex<double>* EvaluatorUtils::randomComplexArray(long size, double bound) {
	complex<double>* res = new complex<double>[size];  // 分配复数数组内存
	for (long i = 0; i < size; ++i) {
		res[i] = randomComplex(bound);                  // 逐个生成随机复数
	}
	return res;
}

// 生成单位圆上的随机复数数组
complex<double>* EvaluatorUtils::randomCircleArray(long size, double bound) {
	complex<double>* res = new complex<double>[size];  // 分配数组内存
	for (long i = 0; i < size; ++i) {
		res[i] = randomCircle(bound);                     // 逐个生成单位圆上的复数
	}
	return res;
}

//----------------------------------------------------------------------------------
//   ROTATIONS  数组旋转实现
//----------------------------------------------------------------------------------

// 数组原地左旋 rotSize 位
// 使用 GCD 分解法：将数组分解为 gcd(rotSize, size) 个独立循环
// 每个循环长度 = size / gcd，循环内逐元素前移
// 此算法时间复杂度 O(size)，空间复杂度 O(1)
void EvaluatorUtils::leftRotateAndEqual(complex<double>* vals, const long size, const long rotSize) {
	long remrotSize = rotSize % size;              // 实际旋转位数（对 size 取模）
	if(remrotSize != 0) {
		long divisor = gcd(remrotSize, size);      // 计算 GCD，决定循环个数
		long steps = size / divisor;               // 每个循环的长度
		for (long i = 0; i < divisor; ++i) {       // 对每个独立循环
			complex<double> tmp = vals[i];         // 保存循环起点元素
			long idx = i;
			for (long j = 0; j < steps - 1; ++j) {  // 循环内逐元素前移
				vals[idx] = vals[(idx + remrotSize) % size];  // 后面的元素前移
				idx = (idx + remrotSize) % size;               // 更新索引
			}
			vals[idx] = tmp;                       // 将起点元素放到循环末尾
		}
	}
}

// 数组原地右旋 rotSize 位
// 右旋等价于左旋 (size - rotSize) 位
void EvaluatorUtils::rightRotateAndEqual(complex<double>* vals, const long size, const long rotSize) {
	long remrotSize = rotSize % size;                    // 实际右旋位数
	long leftremrotSize = (size - remrotSize) % size;    // 转化为左旋位数
	leftRotateAndEqual(vals, size, leftremrotSize);      // 调用左旋实现
}
