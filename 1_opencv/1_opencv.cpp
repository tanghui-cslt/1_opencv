#include <iostream>
#include "Fourier.h"
#include <opencv2\opencv.hpp>

using namespace std;
using namespace cv;

template <class T> T** Matconvert2Array(Mat image);
Mat array2Mat(double **data, int M, int N);

int main(int argc, char** argv)
{
	Mat img = imread("11.jpg");
	Mat grayImg;
	cvtColor(img, grayImg, CV_BGR2GRAY);
	Mat doubleImg;
	grayImg.convertTo(doubleImg, CV_64FC1, 1 / 255.0);

	imshow("img",doubleImg);
	double ** data = Matconvert2Array<double>(grayImg);

	int M = img.rows;
	// cols(width) of the image
	int N = img.cols;
	//double **imageData =DFT(data, M, N);
	double **imageData = DFT(grayImg, M, N);
	Mat newImage = array2Mat(imageData,M,N);

	//imshow("img", newImage);
	//imshow('dimg', doubleImg);
	//imshow('nowe', newImage);
	imshow("newImg1", newImage);
	waitKey(0);

	delete[] data;
	delete[] imageData;
	return 0;
}

template <class T> T **Matconvert2Array(Mat image)
{
	int iRows = image.rows;
	// cols(width) of the image
	int iCols = image.cols;
	int channels = image.channels();
	T ** data = new T*[iRows];
	cout << iRows << " " << iCols <<" "<<channels<< endl;
	//double **ucharData = new double*[iRows];
	
	for (int i = 0; i < iRows; i++)
	{
		//ucharData[i] = new double[iCols];
		data[i] = new T[iCols];
		data[i] = image.ptr<T>(i);
		for (int j = 0; j < iCols; j++)
		{
			if (data[i][j] < 1e-10)
			{
				data[i][j] = 0.0;
			}
			//cout<< data[i][j]<<" ";

		}
		cout <<data[1][i]<< "  ";
	}
	return data;
}

Mat array2Mat(double **data, int M, int N)
{
	Mat image(M, N, CV_64FC1, (double*)data);
	return image;
}