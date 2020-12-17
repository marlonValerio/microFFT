microFFT
==========

Fast Fourier Transform for microcontrollers

This is a fork from https://github.com/kosme/arduinoFFT. That is an Arduino library written in C++ by Enrique Condes. I just adapted it for C projects and give it my own programming style :D.

<del>This is a microcontroller C library for computing FFT.</del> It was tested on ARM microcontrollers.

### Installation

Clone it and include to your project.

`#include "microFTT.h"`

### API

* **FFT_Init**(float *vReal, float *vImag, uint16_t samples, float samplingFrequency);
FFT Initialization
vReal: Is the float array where the signal input samples are stored to be proccesed and, after proccesing, load the real part of the FFT.
vImag: All their values shuld be zero before compute the FFT. Load the imaginary part of the FFT.
sample: Is the number of samples taken in count to the FFT. It should be a multiple of 2.
samplingFrequency: Is the frequency which was used to sampling the signal. It should be 2 Times the maximum signal frequency.
* **FFT_Windowing**(uint8_t windowType, uint8_t dir);
Performs a windowing function on the values array. The possible windowing options are:
    * FFT_WIN_TYP_RECTANGLE
    * FFT_WIN_TYP_HAMMING
    * FFT_WIN_TYP_HANN
    * FFT_WIN_TYP_TRIANGLE
    * FFT_WIN_TYP_NUTTALL
    * FFT_WIN_TYP_BLACKMAN
    * FFT_WIN_TYP_BLACKMAN_NUTTALL
    * FFT_WIN_TYP_BLACKMAN_HARRIS
    * FFT_WIN_TYP_FLT_TOP
    * FFT_WIN_TYP_WELCH
* **Compute**(uint8_t dir);
Calcuates the Fast Fourier Transform.
* **ComplexToMagnitude**(void);
Obtains the magnitude of the frequency from the FFT output.
* **MajorPeak**();
Looks for and returns the frequency of the biggest spike in the analyzed signal.
* **Revision**(void);
Returns the library revision.
