/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAAN_TIMEUTILS_H_  // 头文件保护宏
#define HEAAN_TIMEUTILS_H_

#include "Common.h"  // 引入公共头文件（包含 sys/time.h 提供 timeval）

struct timeval;  // POSIX 时间结构体前置声明（实际定义在 sys/time.h 中）

using namespace std;

// ============================================================================
// 计时工具类
// 使用 POSIX gettimeofday() 实现毫秒级计时
// 在测试中用于测量各操作（加密、解密、同态运算等）的耗时
// ============================================================================
class TimeUtils {
public:

	struct timeval startTime, stopTime;  // 开始/结束时间戳
	double timeElapsed;                  // 记录的耗时（毫秒）

	//-----------------------------------------

	TimeUtils();  // 构造函数

	//-----------------------------------------

	/**
	 * starts timer
	 * @param[in] string message  打印的提示消息（如 "Encrypt batch"）
	 */
	void start(string msg);

	/**
	 * stops timer and prints time elapsed in console
	 * @param[in] string message  打印的提示消息
	 */
	void stop(string msg);

	//-----------------------------------------
};

#endif
