#pragma once
#ifndef _FOURIER_H_
#define _FOURIER_H_
#include <opencv2\opencv.hpp>

double ** DFT(double **data, int M, int N);
double **DFT(cv::Mat &data, int M, int N);

#endif
