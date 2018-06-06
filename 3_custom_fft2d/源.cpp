#include<iostream>
#include <string>
#include<fstream>
#include<algorithm>
#include<cmath>
#include<ctime>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;
//using namespace std;
const double eps(1e-8);

typedef long long lint;

const double PI = acos(-1.0);

class SComplex
{
public:
	double real, image;
	SComplex(double _real, double _image)
	{
		real = _real;
		image = _image;
	}
	SComplex(double value)
	{
		real = value;
		image = value;
	}
	SComplex() {}
	double norm()
	{
		return sqrt(real*real + image * image);
	}
};

SComplex operator + (const SComplex &c1, const SComplex &c2)
{
	return SComplex(c1.real + c2.real, c1.image + c2.image);
}

SComplex operator - (const SComplex &c1, const SComplex &c2)
{
	return SComplex(c1.real - c2.real, c1.image - c2.image);
}

SComplex operator * (const SComplex &c1, const SComplex &c2)
{
	return SComplex(c1.real*c2.real - c1.image*c2.image, c1.real*c2.image + c1.image*c2.real);
}
SComplex operator / (const SComplex &c1, const double data)
{
	return SComplex(c1.real / data, c1.image / data);
}
ostream & operator << (ostream & os, const SComplex &c)
{
	os << "(" << c.real << " " << c.image << ")";
	return os;
}

int getBitNum(int len, int scale)
{
	int i = -1;
	while (len != 0)
	{

		len /= scale;
		//cout<<"len "<<len<< ' ';
		i++;
	}
	return i;
}
int * transBit(int len, int scale)
{
	int *data = new int[len];
	int *new_data = new int[len];

	int times = getBitNum(len, scale);
	int *index_scale = new int[scale];//每个点的新位置
	for (int i = 0; i < scale; i++)
	{
		data[i] = i;
		index_scale[i] = 0;
		new_data[i] = 0;
	}
	for (int i = scale; i < len; i++)
	{
		data[i] = i;
		new_data[i] = 0;
	}

	for (int i = 0; i < times - 1; i++)
	{
		int block = pow(scale, i);//总共有block块
		int interval = len / block;//每一块有interval个元素

		for (int j = 0; j < len; j++)
		{
			int id = j %scale;
			if (j%interval != 0)
			{
				new_data[index_scale[id]] = data[j];//交换数据到对应的位置
				index_scale[id] ++;
				//cout<<j<<" ";
			}
			else
			{
				int seg_num = interval / scale;
				//cout<<endl;
				for (int k = 0; k < scale; k++)
				{
					index_scale[k] = k*seg_num + j; //下标的新位置 ：每一段的位置+起始位置
													//cout<<"i="<< k << " "<< index_scale[k]<<" ";
				}
				new_data[index_scale[id]] = data[j];//交换数据到对应的位置
				index_scale[id] ++;
				//cout<<endl;
			}
			//cout<<j<<" ";
		}

		for (int j = 0; j < len; j++)
		{
			data[j] = new_data[j];
			//cout<<data[j]<<" ";
		}

	}
	delete[] new_data;
	return data;
}
int *rev(int row)
{
	int *ret = new int[row];
	for (int i = 0; i < row; i++)
	{
		ret[0] = 0;
	}
	for (int k = 0; k < row; k++)
	{
		for (int i = 0; (1 << i) < row; i++)
		{
			ret[k] <<= 1;
			if (k & (1 << i)) ret[k] |= 1;
		}
	}
	return ret;
}
// 计算FFT flag = 1表示列 flag = 0 表示行
SComplex ** FFT(SComplex **a, int id, int row, int col, int flag)
{
	SComplex * temp = NULL;
	int len = 0;
	if (flag == 1)
	{
		len = row;
		int *new_index = rev(row);
		temp = new SComplex[row];
		for (int i = 0; i < row; i++)
		{
			//temp[new_index[i]] = a[i][id];
			temp[i] = a[i][id];
		}
	}
	else if (flag == 0)
	{
		len = col;
		temp = new SComplex[col];
		for (int i = 0; i < col; i++)
		{
			temp[i] = a[id][i];
		}

	}

	else
	{
		cout << "error! can't indicate row or col\n";
		exit(1);
	}

	for (int s = 1; (1 << s) <= len; s++)
	{
		int m = (1 << s);
		int interval = (m >> 1);
		SComplex wm = SComplex(cos(-2 * PI / m), sin(-2 * PI / m));
		for (int k = 0; k < len; k += m)
		{
			SComplex w = SComplex(1, 0);
			for (int j = 0; j < interval; j++)
			{
				SComplex t = w*temp[k + j + interval];
				SComplex u = temp[k + j];

				temp[k + j] = u + t;
				temp[k + j + interval] = u - t;
				w = w*wm;
			}
		}
	}

	if (flag == 1)
	{
		for (int i = 0; i < row; i++)
		{
			a[i][id] = temp[i];
		}
	}
	else
	{
		for (int i = 0; i < col; i++)
		{
			a[id][i] = temp[i];
		}
	}

	delete[] temp;
	return a;
}
// 计算FFT flag = 1表示列 flag = 0 表示行
SComplex** FFT_any_base(SComplex** a, int id, int row, int col, int flag, int scale)
{
	SComplex * A = NULL;
	int len = 0;

	if (flag == 1)
	{
		len = row;
		int *new_index = transBit(len, scale);
		A = new SComplex[row];
		for (int i = 0; i < row; i++)
		{
			A[new_index[i]] = a[i][id];
		}
	}
	else if (flag == 0)
	{
		len = col;
		int *new_index = transBit(len, scale);
		A = new SComplex[col];
		for (int i = 0; i < col; i++)
		{
			A[new_index[i]] = a[id][i];
		}

	}

	else
	{
		cout << "error! can't indicate row or col\n";
		exit(1);
	}
	SComplex *t = new SComplex[scale];
	for (int s = 1; pow(scale, s) <= len; s++)
	{
		int m = pow(scale, s);

		int block = len / m;
		//总共有多少块

		for (int k = 0; k < len; k += m)//这一层结点的包含数组元素个数都是(1 << s)
		{

			int temp_value = m / scale;
			double coeff = 0;
			int iter_times = 0;
			//每一块，有多少组

			for (int j = 0; j < (temp_value); j++)//折半引理, 根据两个子节点计算父亲节点
			{
				//wt表示外部的运算。
				//coeff = -2*PI*j*iter_times/len;
				//SComplex wt = SComplex(cos(coeff),sin(coeff));

				iter_times++;
				for (int i = 0; i < scale; i++)
				{
					coeff = -2 * PI*j*i*block / len;
					SComplex wt = SComplex(cos(coeff), sin(coeff));
					t[i] = wt*A[k + j + i*temp_value];

				}
				coeff = coeff*coeff;
				for (int i = 0; i < scale; i++)
				{
					A[k + j + i*temp_value] = SComplex(0, 0);
					for (int l = 0; l < scale; l++)
					{
						double rotate = (-2 * PI*i*l) / scale;
						SComplex coe(cos(rotate), sin(rotate));
						A[k + j + i*temp_value] = A[k + j + i*temp_value] + t[l] * coe;
					}
				}

			}
		}

	}

	if (flag == 1)
	{
		for (int i = 0; i < row; i++)
		{
			a[i][id] = A[i];
		}
	}
	else
	{
		for (int i = 0; i < col; i++)
		{
			a[id][i] = A[i];
		}
	}
	delete[] t;
	delete[] A;
	return a;
}

SComplex** FFT2D(SComplex** a, int row, int col, int &newRowLen, int &newColLen, int scale = 2, int ways = 1)//对长度为len(2的幂)的数组进行DFT变换
{
	double start, end, cost;
	if (ways == 0)	scale = 2;
	//getBitNum(len,scale);
	int BitRowNum = getBitNum(row, scale);
	int BitColNum = getBitNum(col, scale);
	//cout << "Bit = " << BitRowNum << " Bit = " << BitColNum << endl;
	SComplex **A = NULL;
	//  对齐到2的次幂
	newRowLen = pow(scale, BitRowNum);
	newColLen = pow(scale, BitColNum);
	if (newRowLen < row)	BitRowNum += 1;
	if (newColLen < col)	BitColNum += 1;
	newRowLen = pow(scale, BitRowNum);
	newColLen = pow(scale, BitColNum);
	//cout << row << " " << newRowLen << " " << col << " " << newColLen << endl;

	//初始化为0
	A = new SComplex*[newRowLen];
	SComplex c(0, 0);
	for (int i = 0; i < newRowLen; i++)
	{
		A[i] = new SComplex[newColLen];
		for (int j = 0; j < newColLen; j++)
		{
			A[i][j] = c;
		}
	}
	//赋值数据
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			A[i][j] = a[i][j];
		}
	}


	start = clock();
	//cout<<"--------------\n";
	// 标记交换过的行和列，以免交换两次。
	bool *testRow = new bool[newRowLen];
	bool *testCol = new bool[newColLen];
	for (int i = 0; i < newRowLen; i++)
	{
		testRow[i] = true;
	}
	for (int i = 0; i < newColLen; i++)
	{
		testCol[i] = true;
	}
	//计算每一列的fft 
	for (int i = 0; i < newColLen; i++)
	{
		if (ways == 1)
			A = FFT_any_base(A, i, newRowLen, newColLen, 1, scale);
		else
		{
			A = FFT(A, i, newRowLen, newColLen, 1);
			//cout << "----------------\n";
			//cout << i<<" "<<newRowLen << endl;
		}
	}


	//计算每一行的fft
	for (int i = 0; i < newRowLen; i++)
	{
		if (ways == 1)
			A = FFT_any_base(A, i, newRowLen, newColLen, 0, scale);
		else
			A = FFT(A, i, newRowLen, newColLen, 0);
	}
	end = clock();
	cost = (end - start) / CLOCKS_PER_SEC;
	if (ways == 1)
	{
		cout << "-----快速傅里叶变换-------\n" << "row = " << newRowLen << " col = " << newColLen << "  基底= " << scale << " time = " << cost << endl;
	}
	else
	{

		cout << "-----快速傅里叶变换-------\n" << "row = " << newRowLen << " col = " << newColLen << "  基底= " << 2 << " time = " << cost << endl;
	}
	//cout << "快速傅里叶变换时间 = " << cost << endl;
	return A;
}

SComplex * dft1d(SComplex *a, int len, int DFT)
{
	SComplex *A = new SComplex[len];
	SComplex complex_i(0, 1);
	for (int i = 0; i < len; i++)
	{
		A[i] = SComplex(0, 0);
	}
	//cout<<endl;
	for (int i = 0; i < len; i++)
	{
		for (int j = 0; j < len; j++)
		{
			//cout<<(2*PI*i*j/(double)len)<<" ";
			double theta = PI*(-2.0)*(DFT)*(double)(i*j) / (double)len;
			SComplex temp(cos(theta), sin(theta));
			A[i] = A[i] + a[j] * temp;
		}
		//cout<<" "<<i<<endl;
	}
	if (DFT == -1)
	{
		for (int i = 0; i < len; i++)
		{
			A[i] = A[i] / len;
		}
	}

	return A;
}
SComplex ** dft2d_2(SComplex **a, int row, int col, int DFT = 1)
{
	double  start, end, cost;
	start = clock();
	SComplex **A = new SComplex*[row];


	for (int i = 0; i < row; i++)
	{
		//cout<<i<<endl;
		A[i] = dft1d(a[i], col, DFT);
		//cout<<"------------\n";

	}

	SComplex  *temp = new SComplex[row];

	for (int j = 0; j < col; j++)
	{

		SComplex  *fftTemp = NULL;
		for (int i = 0; i < row; i++)
			temp[i] = A[i][j];

		fftTemp = dft1d(temp, row, DFT);

		for (int i = 0; i < row; i++)
			A[i][j] = fftTemp[i];

	}

	end = clock();
	cost = (end - start)*1.0 / CLOCKS_PER_SEC;
	cout << "-----离散傅里叶变换-------\n" << "row = " << row << " col = " << col << " time = " << cost << endl;
	delete[] temp;
	return A;
}
SComplex **Matconvert2Array(Mat image)
{
	int iRows = image.rows;
	int iCols = image.cols;
	int channels = image.channels();
	SComplex ** data = new SComplex*[iRows];
	for (int i = 0; i < iRows; i++)
	{
		double *temp = new double[iCols];
		data[i] = new SComplex[iCols];
		temp = image.ptr<double >(i);
		for (int j = 0; j < iCols; j++)
		{
			data[i][j] = SComplex(temp[j], 0);


		}

	}
	return data;
}

double *** complex2double(SComplex **data, int row, int col)
{
	double *** temp = new double **[2];

	temp[0] = new double *[row];
	temp[1] = new double *[row];
	for (int i = 0; i < row; i++)
	{
		temp[0][i] = new double[col];
		temp[1][i] = new double[col];
		for (int j = 0; j < col; j++)
		{
			temp[0][i][j] = data[i][j].real;
			temp[1][i][j] = data[i][j].image;
		}
	}

	return temp;
}
Mat double2Mat(double **data, int M, int N)
{
	Mat image(M, N, CV_64FC1);
	for (int i = 0; i < M; i++)
	{
		double *temp = image.ptr<double>(i);
		for (int j = 0; j < N; j++)
		{
			temp[j] = data[i][j];
		}
	}


	return image;
}
Mat complex2Mat(SComplex **data, int row, int col)
{

	Mat matData;
	double ***doubleData = complex2double(data, row, col);

	Mat mat1 = double2Mat(doubleData[0], row, col);
	double *temp = mat1.ptr<double>(0);

	Mat mat2 = double2Mat(doubleData[1], row, col);
	Mat planes[] = { Mat_<double>(mat1), Mat_<double>(mat2) };
	merge(planes, 2, matData);

	delete[] doubleData;
	return matData;
}
void display_frequency(SComplex **fft, int row, int col, String name)
{
	Mat complexI = complex2Mat(fft, row, col);
	Mat matData;
	double ***doubleData = complex2double(fft, row, col);

	Mat mat1 = double2Mat(doubleData[0], row, col);
	double *temp = mat1.ptr<double>(0);

	Mat mat2 = double2Mat(doubleData[1], row, col);
	magnitude(mat1, mat2, mat1);
	Mat magI = mat1;
	magI += Scalar::all(1);
	log(magI, magI);                //转换到对数尺度(logarithmic scale)

									//如果有奇数行或列，则对频谱进行裁剪
	magI = magI(Rect(0, 0, magI.cols&-2, magI.rows&-2));

	//重新排列傅里叶图像中的象限，使得原点位于图像中心
	int cx = magI.cols / 2;
	int cy = magI.rows / 2;

	Mat q0(magI, Rect(0, 0, cx, cy));       //左上角图像划定ROI区域
	Mat q1(magI, Rect(cx, 0, cx, cy));      //右上角图像
	Mat q2(magI, Rect(0, cy, cx, cy));      //左下角图像
	Mat q3(magI, Rect(cx, cy, cx, cy));     //右下角图像

											//变换左上角和右下角象限
	Mat tmp;
	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	//变换右上角和左下角象限
	q1.copyTo(tmp);
	q2.copyTo(q1);
	tmp.copyTo(q2);

	//归一化处理，用0-1之间的浮点数将矩阵变换为可视的图像格式
	normalize(magI, magI, 0, 1, CV_MINMAX);

	//outfile << magI << endl;
	imshow(name + "频谱图", magI);
}
int main()
{
	Mat I = imread("14.jpg", IMREAD_GRAYSCALE);
	//ways = 0 专用的fft 基底要设为2
	//ways = 1, 任意基底，基底任意设
	int ways = 1;
	const int row = I.rows, col = I.cols;
	int scale = 4;
	if (I.empty())
	{
		cout << "图像加载失败!" << endl;
		return -1;
	}
	else
		cout << "图像加载成功!" << endl;
	Mat grayImg;
	int channels = I.channels();
	if (channels > 1)
		cvtColor(I, grayImg, CV_BGR2GRAY);
	else
		grayImg = I.clone();
	Mat doubleImg;
	grayImg.convertTo(doubleImg, CV_64FC1, 1);

	SComplex** a = NULL;
	a = Matconvert2Array(doubleImg);
	SComplex **fft = NULL;
	SComplex **dft = NULL;

	int newRow = row, newCol = col;
	//
	dft = dft2d_2(a, row, col, 1);
	display_frequency(dft, row, col, "dft");

	fft = FFT2D(a, row, col, newRow, newCol, scale, ways);
	display_frequency(fft, newRow, newCol, "fft");

	imshow("输入图像", I);

	waitKey(0);
	I.release();
	grayImg.release();
	doubleImg.release();
	//outfile.close();
	delete[]  a;
	delete[] fft;
	delete[] dft;
	return 0;
}
