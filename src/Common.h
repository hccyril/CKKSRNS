/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_COMMON_H_   // 头文件保护宏，防止重复包含
#define HEAANNTT_COMMON_H_

// 公共头文件：集中引入整个项目所需的标准库
// 所有其他源文件通过 #include "Common.h" 获得基础类型支持

#include <iostream>           // 标准输入输出流（cout, endl 等），用于打印日志与调试信息
#include <stdio.h>            // C 标准 I/O 库
#include <stdint.h>           // 固定宽度整数类型定义：uint64_t, uint32_t 等
                              // 全方案的多项式系数、密文分量均以 uint64_t 存储

#include <stdexcept>          // 异常类：invalid_argument 等，用于参数校验失败时抛出异常
#include <vector>             // 动态数组容器
#include <sys/time.h>         // POSIX 计时函数 gettimeofday()，TimeUtils 类使用
#include <string>             // 字符串类，用于函数名标识、日志输出等
#include <math.h>             // 数学函数库：sqrt, log, floor, ceil, cos, sin 等
                              // 编码/解码中的 FFT、高斯采样等需要

#endif
