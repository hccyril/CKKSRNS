/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "StringUtils.h"  // 引入字符串工具类头文件

//----------------------------------------------------------------------------------
//   SHOW ARRAY  数组打印实现
//----------------------------------------------------------------------------------

// 打印 uint64_t 数组，格式: [val0, val1, ..., valN]
void StringUtils::show(uint64_t* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";  // 逐元素打印，逗号分隔
	}
	cout << "]" << endl;
}

// 打印 long 数组
void StringUtils::show(long* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";
	}
	cout << "]" << endl;
}

// 打印 double 数组
void StringUtils::show(double* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";
	}
	cout << "]" << endl;
}

// 打印复数数组
void StringUtils::show(complex<double>* vals, long size) {
	cout << "[";
	for (long i = 0; i < size; ++i) {
		cout << vals[i] << ", ";
	}
	cout << "]" << endl;
}


//----------------------------------------------------------------------------------
//   SHOW & COMPARE ARRAY  打印并对比
//----------------------------------------------------------------------------------


// 对比两个 double 标量：打印明文值、解密值、误差
void StringUtils::showcompare(double val1, double val2, string prefix) {
	cout << "---------------------" << endl;
	cout << "m" + prefix + ":" << val1 << endl;     // m = message（明文值）
	cout << "d" + prefix + ":" << val2 << endl;     // d = decrypted（解密值）
	cout << "e" + prefix + ":" << val1-val2 << endl; // e = error（误差 = 明文 - 解密）
	cout << "---------------------" << endl;
}

// 对比两个复数标量
void StringUtils::showcompare(complex<double> val1, complex<double> val2, string prefix) {
	cout << "---------------------" << endl;
	cout << "m" + prefix + ":" << val1 << endl;
	cout << "d" + prefix + ":" << val2 << endl;
	cout << "e" + prefix + ":" << val1-val2 << endl;
	cout << "---------------------" << endl;
}

// 逐元素对比两个 double 数组
void StringUtils::showcompare(double* vals1, double* vals2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;    // 第 i 个明文值
		cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;    // 第 i 个解密值
		cout << "e" + prefix + ": " << i << " :" << (vals1[i]-vals2[i]) << endl;  // 第 i 个误差
		cout << "---------------------" << endl;
	}
}

// 逐元素对比两个复数数组
void StringUtils::showcompare(complex<double>* vals1, complex<double>* vals2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
		cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
		cout << "e" + prefix + ": " << i << " :" << (vals1[i]-vals2[i]) << endl;
		cout << "---------------------" << endl;
	}
}


// double 数组 vs double 标量
void StringUtils::showcompare(double* vals1, double val2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
		cout << "d" + prefix + ": " << i << " :" << val2 << endl;
		cout << "e" + prefix + ": " << i << " :" << (vals1[i]-val2) << endl;
		cout << "---------------------" << endl;
	}
}

// 复数数组 vs 复数标量
void StringUtils::showcompare(complex<double>* vals1, complex<double> val2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << vals1[i] << endl;
		cout << "d" + prefix + ": " << i << " :" << val2 << endl;
		cout << "e" + prefix + ": " << i << " :" << (vals1[i]-val2) << endl;
		cout << "---------------------" << endl;
	}
}

// double 标量 vs double 数组
void StringUtils::showcompare(double val1, double* vals2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << val1 << endl;
		cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
		cout << "e" + prefix + ": " << i << " :" << (val1-vals2[i]) << endl;
		cout << "---------------------" << endl;
	}
}

// 复数标量 vs 复数数组
void StringUtils::showcompare(complex<double> val1, complex<double>* vals2, long size, string prefix) {
	for (long i = 0; i < size; ++i) {
		cout << "---------------------" << endl;
		cout << "m" + prefix + ": " << i << " :" << val1 << endl;
		cout << "d" + prefix + ": " << i << " :" << vals2[i] << endl;
		cout << "e" + prefix + ": " << i << " :" << (val1-vals2[i]) << endl;
		cout << "---------------------" << endl;
	}
}
