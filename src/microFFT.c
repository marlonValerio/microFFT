/*

	FFT libray
	Copyright (C) 2010 Didier Longueville
	Copyright (C) 2014 Enrique Condes

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "microFFT.h"

/*Mathematial constants*/
#define twoPi 	6.28318531
#define fourPi 	12.56637061
#define sixPi 	18.84955593

#define FFT_LIB_REV 0x10

/* Private Variables */
static uint16_t _samples;
static float 	_samplingFrequency;
static float*	_vReal;
static float*	_vImag;
static uint8_t 	_power;

/* Private Functions and Macros*/
static void 	Swap(float *x, float *y);
static uint8_t 	Exponent(uint16_t value);
#define sq(x) 	((x)*(x))

void FFT_Init(float *vReal, float *vImag, uint16_t samples, float samplingFrequency)
{
	_vReal = vReal;
	_vImag = vImag;
	_samples = samples;
	_samplingFrequency = samplingFrequency;
	_power = Exponent(samples);
}

void FFT_Compute(uint8_t dir)
{// Computes in-place complex-to-complex FFT /
	// Reverse bits /
	uint16_t j = 0;
	for (uint16_t i = 0; i < (_samples - 1); i++) {
		if (i < j) {
			Swap(&_vReal[i], &_vReal[j]);
			if(dir==FFT_REVERSE)
				Swap(&_vImag[i], &_vImag[j]);
		}
		uint16_t k = (_samples >> 1);
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}
	// Compute the FFT  /
	float c1 = -1.0;
	float c2 = 0.0;
	uint16_t l2 = 1;
	for (uint8_t l = 0; (l < _power); l++) {
		uint16_t l1 = l2;
		l2 <<= 1;
		float u1 = 1.0;
		float u2 = 0.0;
		for (j = 0; j < l1; j++) {
			 for (uint16_t i = j; i < _samples; i += l2) {
					uint16_t i1 = i + l1;
					float t1 = u1 * _vReal[i1] - u2 * _vImag[i1];
					float t2 = u1 * _vImag[i1] + u2 * _vReal[i1];
					_vReal[i1] = _vReal[i] - t1;
					_vImag[i1] = _vImag[i] - t2;
					_vReal[i] += t1;
					_vImag[i] += t2;
			 }
			 float z = ((u1 * c1) - (u2 * c2));
			 u2 = ((u1 * c2) + (u2 * c1));
			 u1 = z;
		}
		c2 = sqrt((1.0 - c1) / 2.0);
		c1 = sqrt((1.0 + c1) / 2.0);
		if (dir == FFT_FORWARD) {
			c2 = -c2;
		}
	}
	// Scaling for reverse transform /
	if (dir != FFT_FORWARD) {
		for (uint16_t i = 0; i < _samples; i++) {
			 _vReal[i] /= _samples;
			 _vImag[i] /= _samples;
		}
	}
}

void FFT_ComplexToMagnitude(void)
{ // vM is half the size of vReal and vImag
	for (uint16_t i = 0; i < _samples; i++) {
		_vReal[i] = sqrt(sq(_vReal[i]) + sq(_vImag[i]));
	}
}

void FFT_Windowing(uint8_t windowType, uint8_t dir)
{// Weighing factors are computed once before multiple use of FFT
// The weighing function is symetric; half the weighs are recorded
	float samplesMinusOne = ((float)(_samples) - 1.0);
	for (uint16_t i = 0; i < (_samples >> 1); i++) {
		float indexMinusOne = (float)(i);
		float ratio = (indexMinusOne / samplesMinusOne);
		float weighingFactor = 1.0;
		// Compute and record weighting factor
		switch (windowType) {
		case FFT_WIN_TYP_RECTANGLE: // rectangle (box car)
			weighingFactor = 1.0;
			break;
		case FFT_WIN_TYP_HAMMING: // hamming
			weighingFactor = 0.54 - (0.46 * cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_HANN: // hann
			weighingFactor = 0.54 * (1.0 - cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_TRIANGLE: // triangle (Bartlett)
			weighingFactor = 1.0 - ((2.0 * abs(indexMinusOne - (samplesMinusOne / 2.0))) / samplesMinusOne);
			break;
		case FFT_WIN_TYP_NUTTALL: // nuttall
			weighingFactor = 0.355768 - (0.487396 * (cos(twoPi * ratio))) + (0.144232 * (cos(fourPi * ratio))) - (0.012604 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN: // blackman
			weighingFactor = 0.42323 - (0.49755 * (cos(twoPi * ratio))) + (0.07922 * (cos(fourPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_NUTTALL: // blackman nuttall
			weighingFactor = 0.3635819 - (0.4891775 * (cos(twoPi * ratio))) + (0.1365995 * (cos(fourPi * ratio))) - (0.0106411 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_HARRIS: // blackman harris
			weighingFactor = 0.35875 - (0.48829 * (cos(twoPi * ratio))) + (0.14128 * (cos(fourPi * ratio))) - (0.01168 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_FLT_TOP: // flat top
			weighingFactor = 0.2810639 - (0.5208972 * cos(twoPi * ratio)) + (0.1980399 * cos(fourPi * ratio));
			break;
		case FFT_WIN_TYP_WELCH: // welch
			weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) / (samplesMinusOne / 2.0));
			break;
		}
		if (dir == FFT_FORWARD) {
			_vReal[i] *= weighingFactor;
			_vReal[_samples - (i + 1)] *= weighingFactor;
		}
		else {
			_vReal[i] /= weighingFactor;
			_vReal[_samples - (i + 1)] /= weighingFactor;
		}
	}
}

void FFT_MajorPeak(float* mag_out, float* freq_out, float magFact)
{
	float maxY = 0;
	uint16_t IndexOfMaxY = 0;

	for (uint16_t i = 1; i < ((_samples >> 1) + 1); i++) {
		if ((_vReal[i-1] < _vReal[i]) && (_vReal[i] > _vReal[i+1])) {
			if (_vReal[i] > maxY) {
				maxY = _vReal[i];
				IndexOfMaxY = i;
			}
		}
	}
	float delta = 0.5 * ((_vReal[IndexOfMaxY-1] - _vReal[IndexOfMaxY+1]) / (_vReal[IndexOfMaxY-1] - (2.0 * _vReal[IndexOfMaxY]) + _vReal[IndexOfMaxY+1]));
	float interpolatedX = ((IndexOfMaxY + delta)  * _samplingFrequency) / (_samples-1);
	if(IndexOfMaxY==(_samples >> 1))
		interpolatedX = ((IndexOfMaxY + delta)  * _samplingFrequency) / (_samples);

	*mag_out = _vReal[IndexOfMaxY]/magFact;
	*freq_out = interpolatedX;
}

uint8_t FFT_LibRevision(void)
{
	return(FFT_LIB_REV);
}

uint8_t Exponent(uint16_t value)
{
	uint8_t result = 0;
	while (((value >> result) & 1) != 1) result++;
	return(result);
}

void Swap(float *x, float *y)
{
	float temp = *x;
	*x = *y;
	*y = temp;
}
