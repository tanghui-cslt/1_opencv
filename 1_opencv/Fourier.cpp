#include "Fourier.h"
#include <iostream>
#include <complex>
#include <cmath>

using namespace std;
using namespace  cv;
#define PI 3.14159265358
// 计算指数
// u ,v 表示周期, x,y 表示图像坐标, M,N 表示图像尺寸
complex<double > exp_w(int u, int v, int x, int y, int M, int N)
{
	double theta = 2.0*PI*(1.0*u*x / M + 1.0 *v *y / N);
	return complex<double>(cos(theta), -sin(theta));
}
double **DFT(Mat &data, int M, int N)
{
	complex<double > **temp = new complex<double> *[M];//复数
	double **uniform_double = new double *[M];//保存对数化之后的数据
	double **newData = new double *[M];//图像数据
	double minValue = 255;
	double maxValue = 0;
	for (int i = 0; i < M; i++)
	{
		temp[i] = new complex<double>[N];
		uniform_double[i] = new double[N];
		newData[i] = new double[N];
	}
	for (int i = 0; i < M; i++)
	{
		//cout << i << endl;
		for (int j = 0; j < N; j++)
		{
			temp[i][j] = (0, 0);// init complex 
			for (int l = 0; l < M; l++)
			{
				for (int k = 0; k < N; k++)
				{
					temp[i][j] += data.at<uchar>(l,k)*1.0 * exp_w(i, j, l, k, M, N);
				}
			}
			//cout << temp[i][j] << " ";
			double norm_amplitude = log(sqrt(norm(temp[i][j])) + 1);
			minValue = minValue > norm_amplitude ? norm_amplitude : minValue;
			maxValue = maxValue < norm_amplitude ? norm_amplitude : maxValue;
			//cout << norm_amplitude << " ";
			uniform_double[i][j] = norm_amplitude;
		}
	}
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double temp_255 = (uniform_double[i][j] - minValue) / (maxValue - minValue);
			newData[i][j] = temp_255;

			//if (temp_255 >250 || temp_255< 2)
			{
				//cout << (int)newData[i][j] << " ";
			}
		}
		cout << (int)newData[1][i] << " ";
	}
	delete[] temp;
	delete[] uniform_double;
	return newData;
}

double ** DFT(double **data, int M, int N)
{
	
	complex<double > **temp = new complex<double> *[M];//复数
	double **uniform_double = new double *[M];//保存对数化之后的数据
	double **newData = new double *[M];//图像数据
	double minValue = 255;
	double maxValue = 0;
	for (int i = 0; i < M; i++)
	{
		temp[i] = new complex<double>[N];
		uniform_double[i] = new double[N];
		newData[i] = new double[N];
	}

	cout << "--------------complex----------\n";
	for (int i = 0; i < M; i++)
	{
		//cout << i << endl;
		for (int j = 0; j < N; j++)
		{
			temp[i][j] = (0, 0);// init complex 
			for (int l = 0; l < M; l++)
			{
				for (int k = 0; k < N; k++)
				{
					temp[i][j] += data[l][k] * exp_w(i, j, l, k, M, N);
				}
			}
			//cout << temp[i][j] << " ";
			double norm_amplitude = log(sqrt(norm(temp[i][j]))+1);
			//if (norm_amplitude >6)
			//{
			////	cout << "------" << i << " *****  " << j<<" " << sqrt(norm(temp[i][j]))<<" " << log(sqrt(norm(temp[i][j]) )+1) << endl;
			//	//getchar();
			//}
			minValue = minValue > norm_amplitude ? norm_amplitude : minValue;
			maxValue = maxValue < norm_amplitude ? norm_amplitude : maxValue;
			//cout << norm_amplitude << " ";
			uniform_double[i][j] = norm_amplitude;
			//newData[i][j] = (uchar)norm_amplitude;
		}
		//cout << endl;
	}
	cout << maxValue << "------------------ " << minValue << endl;
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			double temp_255 = (uniform_double[i][j] - minValue) /( maxValue - minValue);
			newData[i][j] = temp_255;

			//if (temp_255 >250 || temp_255< 2)
			{
				//cout << (int)newData[i][j] << " ";
			}
		}
		cout << (int )newData[1][i]<< " ";
	}
	delete[] temp;
	delete[] uniform_double;
	return newData;
}


