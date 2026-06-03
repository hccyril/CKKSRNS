/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

// ============================================================================
// SchemeAlgo.cpp —— 高阶算法实现
// 基于 Scheme 的基本同态运算构建复杂算法:
//   - 幂运算: 反复平方 + rescale
//   - 乘积: 树状乘法
//   - 求逆: 1/(1-x) 的截断几何级数
//   - 指数/sigmoid: Taylor 级数展开
//   - 同态 FFT: 在密文数组上执行 FFT
// ============================================================================
#include "SchemeAlgo.h"

//----------------------------------------------------------------------------------
//   ARRAY ENCRYPTION & DECRYPTION
//----------------------------------------------------------------------------------


Ciphertext* SchemeAlgo::encryptSingleArray(complex<double>* vals, long size) {
	Ciphertext* res = new Ciphertext[size];
	for (long i = 0; i < size; ++i) {
		res[i] = scheme.encryptSingle(vals[i], scheme.context.L);
	}
	return res;
}

Ciphertext* SchemeAlgo::encryptSingleArray(double* vals, long size) {
	Ciphertext* res = new Ciphertext[size];
	for (long i = 0; i < size; ++i) {
		res[i] = scheme.encryptSingle(vals[i], scheme.context.L);
	}
	return res;
}

complex<double>* SchemeAlgo::decryptSingleArray(SecretKey& secretKey, Ciphertext* ciphers, long size) {
	complex<double>* res = new complex<double>[size];
	for (long i = 0; i < size; ++i) {
		res[i] = scheme.decryptSingle(secretKey, ciphers[i]);
	}
	return res;
}

//----------------------------------------------------------------------------------
//   POWERS & PRODUCTS
//----------------------------------------------------------------------------------


// ============================================================================
// powerOf2 —— 计算 cipher^{2^logDegree}
// 通过反复平方实现: cipher → cipher² → cipher⁴ → ... → cipher^{2^logDegree}
// 每次平方后执行 rescale 将 scale 从 Δ² 恢复到 Δ
// 总消耗: logDegree 次乘法 + logDegree 次 rescale
// 消耗 logDegree 层 level
// ============================================================================
Ciphertext SchemeAlgo::powerOf2(Ciphertext& cipher, const long logDegree) {
	Ciphertext res = cipher;
	for (long i = 0; i < logDegree; ++i) {
		scheme.squareAndEqual(res);
		scheme.reScaleByAndEqual(res, 1);
	}
	return res;
}

Ciphertext* SchemeAlgo::powerOf2Extended(Ciphertext& cipher, const long logDegree) {
	Ciphertext* res = new Ciphertext[logDegree + 1];
	res[0] = cipher;
	for (long i = 1; i < logDegree + 1; ++i) {
		res[i] = scheme.square(res[i-1]);
		scheme.reScaleByAndEqual(res[i], 1);
	}
	return res;
}

Ciphertext SchemeAlgo::power(Ciphertext& cipher, const long degree) {
	long logDegree = log2((double)degree);
	long po2Degree = 1 << logDegree;

	Ciphertext res = powerOf2(cipher, logDegree);
	long remDegree = degree - po2Degree;
	if(remDegree > 0) {
		Ciphertext tmp = power(cipher, remDegree);
		scheme.modDownToAndEqual(tmp, tmp.l - 1);
		scheme.multAndEqual(res, tmp);
		scheme.reScaleByAndEqual(res, 1);
	}
	return res;
}

Ciphertext* SchemeAlgo::powerExtended(Ciphertext& cipher, const long degree) {
	Ciphertext* res = new Ciphertext[degree];
	long logDegree = log2((double)degree);
	Ciphertext* cpows = powerOf2Extended(cipher, logDegree);
	long idx = 0;
	for (long i = 0; i < logDegree; ++i) {
		long powi = (1 << i);
		res[idx++] = cpows[i];
		for (int j = 0; j < powi-1; ++j) {
			res[idx] = scheme.modDownTo(res[j], cpows[i].l);
			scheme.multAndEqual(res[idx], cpows[i]);
			scheme.reScaleByAndEqual(res[idx++], 1);
		}
	}
	res[idx++] = cpows[logDegree];
	long degree2 = (1 << logDegree);
	for (int i = 0; i < (degree - degree2); ++i) {
		res[idx] = scheme.modDownTo(res[i], cpows[logDegree].l);
		scheme.multAndEqual(res[idx], cpows[logDegree]);
		scheme.reScaleByAndEqual(res[idx++], 1);
	}
	return res;
}

Ciphertext SchemeAlgo::prodOfPo2(Ciphertext* ciphers, const long logDegree) {
	Ciphertext* res = ciphers;
	for (long i = logDegree - 1; i >= 0; --i) {
		long powih = (1 << i);
		Ciphertext* tmp = new Ciphertext[powih];
		for (long j = 0; j < powih; ++j) {
			tmp[j] = scheme.mult(res[2 * j], res[2 * j + 1]);
			scheme.reScaleByAndEqual(tmp[j], 1);
		}
		res = tmp;
	}
	return res[0];
}

Ciphertext SchemeAlgo::prod(Ciphertext* ciphers, const long degree) {
	long logDegree = log2((double)degree) + 1;
	long idx = 0;
	bool isinit = false;
	Ciphertext res;
	for (long i = 0; i < logDegree; ++i) {
		if(degree & 1 << i) {
			long powi = (1 << i);
			Ciphertext* tmp = new Ciphertext[powi];
			for (long j = 0; j < powi; ++j) {
				tmp[j] = ciphers[idx + j];
			}
			Ciphertext iprod = prodOfPo2(tmp, i);
			if(isinit) {
				long dl = res.l - iprod.l;
				scheme.modDownByAndEqual(res, dl);
				scheme.multAndEqual(res, iprod);
				scheme.reScaleByAndEqual(res, 1);
			} else {
				res = iprod;
				isinit = true;
			}
			idx += powi;
		}
	}
	return res;
}


//----------------------------------------------------------------------------------
//   METHODS ON ARRAYS OF CIPHERTEXTS
//----------------------------------------------------------------------------------


Ciphertext SchemeAlgo::sum(Ciphertext* ciphers, const long size) {
	Ciphertext res = ciphers[0];
	for (long i = 1; i < size; ++i) {
		scheme.addAndEqual(res, ciphers[i]);
	}
	return res;
}

Ciphertext SchemeAlgo::distance(Ciphertext& cipher1, Ciphertext& cipher2) {
	Ciphertext cres = scheme.sub(cipher1, cipher2);
	scheme.squareAndEqual(cres);
	scheme.reScaleByAndEqual(cres, 1);
	partialSlotsSumAndEqual(cres, cres.slots);
	return cres;
}

Ciphertext* SchemeAlgo::multVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size) {
	Ciphertext* res = new Ciphertext[size];
	for (long i = 0; i < size; ++i) {
		res[i] = scheme.mult(ciphers1[i], ciphers2[i]);
	}
	return res;
}

void SchemeAlgo::multAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size) {
	for (long i = 0; i < size; ++i) {
		scheme.multAndEqual(ciphers1[i], ciphers2[i]);
	}
}


Ciphertext* SchemeAlgo::multAndModSwitchVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size) {
	Ciphertext* res = new Ciphertext[size];
	for (long i = 0; i < size; ++i) {
		res[i] = scheme.mult(ciphers1[i], ciphers2[i]);
		scheme.reScaleByAndEqual(res[i], 1);
	}
	return res;
}

void SchemeAlgo::multModSwitchAndEqualVec(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size) {
	for (long i = 0; i < size; ++i) {
		scheme.multAndEqual(ciphers1[i], ciphers2[i]);
		scheme.reScaleByAndEqual(ciphers1[i], 1);
	}
}

Ciphertext SchemeAlgo::innerProd(Ciphertext* ciphers1, Ciphertext* ciphers2, const long size) {
	Ciphertext cip = scheme.mult(ciphers1[size-1], ciphers2[size-1]);

	for (long i = 0; i < size - 1; ++i) {
		Ciphertext cprodi = scheme.mult(ciphers1[i], ciphers2[i]);
		scheme.addAndEqual(cip, cprodi);
	}

	scheme.reScaleByAndEqual(cip, 1);
	return cip;
}

Ciphertext SchemeAlgo::partialSlotsSum(Ciphertext& cipher, const long slots) {
	Ciphertext res = cipher;
	for (long i = 1; i < slots; i <<= 1) {
		Ciphertext rot = scheme.leftRotateFast(res, i);
		scheme.addAndEqual(res, rot);
	}
	return res;
}

void SchemeAlgo::partialSlotsSumAndEqual(Ciphertext& cipher, const long slots) {
	for (long i = 1; i < slots; i <<= 1) {
		Ciphertext rot = scheme.leftRotateFast(cipher, i);
		scheme.addAndEqual(cipher, rot);
	}
}


//----------------------------------------------------------------------------------
//   FUNCTIONS
//----------------------------------------------------------------------------------


// ============================================================================
// inverse —— 求逆运算: 1/cipher
// 基于几何级数: 1/(1-x) = 1 + x + x² + x³ + ...
// 算法（Newton 迭代思路）:
//   初始化: cbar = 1 - cipher, res = 1 + (1-x) = 2-x
//   第 i 步: cpow = cpow² (x^{2^i}), res = res * (1 + cpow)
//   结果收敛到 1/(1-x) 的截断级数
//
// 前提: |cipher| < 1（输入需在单位圆内）
// 精度: 2^steps 次项后的截断误差
// ============================================================================
Ciphertext SchemeAlgo::inverse(Ciphertext& cipher, const long steps) {
	Ciphertext cbar = scheme.negate(cipher);
	scheme.addConstAndEqual(cbar, 1.0);
	Ciphertext cpow = cbar;
	Ciphertext tmp = scheme.addConst(cbar, 1.0);
	scheme.modDownByAndEqual(tmp, 1);
	Ciphertext res = tmp;

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, 1);
		tmp = cpow;
		scheme.addConstAndEqual(tmp, 1.0);
		scheme.multAndEqual(tmp, res);
		scheme.reScaleByAndEqual(tmp, 1);
		res = tmp;
	}
	return res;
}

Ciphertext* SchemeAlgo::inverseExtended(Ciphertext& cipher, const long steps) {
	Ciphertext* res = new Ciphertext[steps];
	Ciphertext cpow = cipher;
	Ciphertext tmp = scheme.addConst(cipher, 1.0);
	scheme.modDownByAndEqual(tmp, 1);
	res[0] = tmp;

	for (long i = 1; i < steps; ++i) {
		scheme.squareAndEqual(cpow);
		scheme.reScaleByAndEqual(cpow, 1);
		tmp = cpow;
		scheme.addConstAndEqual(tmp, 1.0);
		scheme.multAndEqual(tmp, res[i - 1]);
		scheme.reScaleByAndEqual(tmp, 1);
		res[i] = tmp;
	}
	return res;
}

// ============================================================================
// exponent —— 指数函数: exp(cipher)
// 使用 Taylor 级数展开: exp(x) = 1 + x + x²/2! + x³/3! + ...
// 步骤:
//   1. powerExtended: 计算 [cipher, cipher², ..., cipher^degree]
//   2. 从 taylorCoeffsMap 取 Taylor 系数
//   3. 线性组合: res = Σ coeffs[i+1] * cipher^i + coeffs[0]
//   4. rescale: 恢复 scale
//
// 论文 §5 函数求值的核心应用之一
// ============================================================================
Ciphertext SchemeAlgo::exponent(Ciphertext& cipher, long degree) {
	Ciphertext* cpows = powerExtended(cipher, degree);

	double* coeffs = scheme.context.taylorCoeffsMap.at(EXPONENT);

	Ciphertext res = scheme.multByConst(cpows[0], coeffs[1]);
	scheme.addP2AndEqual(res); // 2 lvl
	for (int i = 1; i < degree; ++i) {
		if(abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext aixi = scheme.multByConst(cpows[i], coeffs[i + 1]);
			scheme.modDownToAndEqual(res, aixi.l);
			scheme.addAndEqual(res, aixi);
		}
	}
	scheme.reScaleByAndEqual(res, 1);
	return res;
}

// ============================================================================
// sigmoid —— Sigmoid 函数: σ(x) = exp(x)/(1+exp(x)) = 1/(1+e^{-x})
// 使用 Taylor 级数近似（系数预存在 taylorCoeffsMap 中）
// 实现方式与 exponent 相同，仅系数不同
// 论文 §5: sigmoid 在机器学习场景中特别重要
// ============================================================================
Ciphertext SchemeAlgo::sigmoid(Ciphertext& cipher, long degree) {
	Ciphertext* cpows = powerExtended(cipher, degree);

	double* coeffs = scheme.context.taylorCoeffsMap.at(SIGMOID);

	Ciphertext res = scheme.multByConst(cpows[0], coeffs[1]);
	scheme.addP2hAndEqual(res); // 2 lvl
	for (int i = 1; i < degree; ++i) {
		if(abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext aixi = scheme.multByConst(cpows[i], coeffs[i + 1]);
			scheme.modDownToAndEqual(res, aixi.l);
			scheme.addAndEqual(res, aixi);
		}
	}
	scheme.reScaleByAndEqual(res, 1);
	return res;
}

Ciphertext SchemeAlgo::function(Ciphertext& cipher, string& funcName, const long degree) {
	Ciphertext* cpows = powerExtended(cipher, degree);

	double* coeffs = scheme.context.taylorCoeffsMap.at(funcName);

	Ciphertext res = scheme.multByConst(cpows[0], coeffs[1]);
	scheme.addConstAndEqual(res, coeffs[0]); // 2 lvl

	for (int i = 1; i < degree; ++i) {
		if(abs(coeffs[i + 1]) > 1e-27) {
			Ciphertext aixi = scheme.multByConst(cpows[i], coeffs[i + 1]);
			scheme.modDownToAndEqual(res, aixi.l);
			scheme.addAndEqual(res, aixi);
		}
	}
	scheme.reScaleByAndEqual(res, 1);
	return res;
}

Ciphertext* SchemeAlgo::functionExtended(Ciphertext& cipher, string& funcName, const long degree) {
	Ciphertext* cpows = powerExtended(cipher, degree);

	double* coeffs = scheme.context.taylorCoeffsMap.at(funcName);

	Ciphertext aixi = scheme.multByConst(cpows[0], coeffs[1]);

//	scheme.addConstAndEqual(aixi, coeffs[0], dlogp);

	Ciphertext* res = new Ciphertext[degree];
	res[0] = aixi;
	for (long i = 1; i < degree; ++i) {
		if(abs(coeffs[i + 1]) > 1e-27) {
			aixi = scheme.multByConst(cpows[i], coeffs[i + 1]);
			Ciphertext ctmp = scheme.modDownTo(res[i - 1], aixi.l);
			scheme.addAndEqual(aixi, ctmp);
			res[i] = aixi;
		} else {
			res[i] = res[i - 1];
		}
	}
	for (long i = 0; i < degree; ++i) {
		scheme.reScaleByAndEqual(res[i], 1);
	}
	return res;
}


//----------------------------------------------------------------------------------
//   FFT & FFT INVERSE
//----------------------------------------------------------------------------------


void SchemeAlgo::bitReverse(Ciphertext* ciphers, const long size) {
	for (long i = 1, j = 0; i < size; ++i) {
		long bit = size >> 1;
		for (; j >= bit; bit>>=1) {
			j -= bit;
		}
		j += bit;
		if(i < j) {
			swap(ciphers[i], ciphers[j]);
		}
	}
}

void SchemeAlgo::fft(Ciphertext* ciphers, const long size) {
	bitReverse(ciphers, size);
	for (long len = 2; len <= size; len <<= 1) {
		long shift = scheme.context.M / len;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < len / 2; ++j) {
				Ciphertext u = ciphers[i + j];
				scheme.multByMonomialAndEqual(ciphers[i + j + len / 2], shift * j);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + len / 2]);
				scheme.sub2AndEqual(u, ciphers[i + j + len / 2]);
			}
		}
	}
}

void SchemeAlgo::fftInvLazy(Ciphertext* ciphers, const long size) {
	bitReverse(ciphers, size);
	for (long len = 2; len <= size; len <<= 1) {
		long shift = scheme.context.M - scheme.context.M / len;
		for (long i = 0; i < size; i += len) {
			for (long j = 0; j < len / 2; ++j) {
				Ciphertext u = ciphers[i + j];
				scheme.multByMonomialAndEqual(ciphers[i + j + len / 2], shift * j);
				scheme.addAndEqual(ciphers[i + j], ciphers[i + j + len / 2]);
				scheme.sub2AndEqual(u, ciphers[i + j + len / 2]);
			}
		}
	}
}

void SchemeAlgo::fftInv(Ciphertext* ciphers, const long size) {
	fftInvLazy(ciphers, size);

}
