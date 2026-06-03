/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "TimeUtils.h"  // 引入计时工具类头文件

// 构造函数：初始化累计耗时为 0
TimeUtils::TimeUtils() {
	timeElapsed = 0;  // 初始化耗时计数器
}

// 开始计时：记录当前时间并打印开始消息
// gettimeofday 获取当前系统时间（秒 + 微秒精度）
void TimeUtils::start(string msg) {
	cout << "------------------" << endl;
	cout <<"Start " + msg << endl;       // 打印操作名称
	gettimeofday(&startTime, 0);         // 记录开始时间
}

// 结束计时：计算耗时并打印
// 耗时 = (stopTime - startTime)，转换为毫秒
void TimeUtils::stop(string msg) {
	gettimeofday(&stopTime, 0);                                              // 记录结束时间
	timeElapsed = (stopTime.tv_sec - startTime.tv_sec) * 1000.0;             // 秒数差 × 1000 → 毫秒
	timeElapsed += (stopTime.tv_usec - startTime.tv_usec) / 1000.0;          // 微秒差 / 1000 → 毫秒
	cout << msg +  " time = "<< timeElapsed << " ms" << endl;                // 打印耗时
	cout << "------------------" << endl;
}
