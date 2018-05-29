#include <iostream>
#include <fstream>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

int main()
{
	Mat I = imread("13.jpg", IMREAD_GRAYSCALE);       //����ͼ��Ҷ�ͼ
	ofstream outfile("data.txt");
	ofstream outfile1("result.txt");
	
														//�ж�ͼ���Ƿ���سɹ�
	if (I.empty())
	{
		cout << "ͼ�����ʧ��!" << endl;
		return -1;
	}
	else
		cout << "ͼ����سɹ�!" << endl << endl;

	Mat padded;                 //��0�������ͼ�����
	int m = getOptimalDFTSize(I.rows);
	int n = getOptimalDFTSize(I.cols);
	//m = 256, n = 256;
	cout << " m = " << m << " n = " << n << endl;
	//�������ͼ��I���������Ϊpadded���Ϸ����󷽲�����䴦��
	copyMakeBorder(I, padded, 0, m - I.rows, 0, n - I.cols, BORDER_CONSTANT, Scalar::all(0));
	outfile << padded << endl;
	outfile.close();
	Mat planes[] = { Mat_<float>(padded), Mat::zeros(padded.size(),CV_32F) };
	Mat complexI;
	merge(planes, 2, complexI);     //��planes�ںϺϲ���һ����ͨ������complexI

	dft(complexI, complexI);        //���и���Ҷ�任
	
									//�����ֵ��ת���������߶�(logarithmic scale)
									//=> log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
	split(complexI, planes);        //planes[0] = Re(DFT(I),planes[1] = Im(DFT(I))
									//��planes[0]Ϊʵ��,planes[1]Ϊ�鲿
	magnitude(planes[0], planes[1], planes[0]);     //planes[0] = magnitude
	Mat magI = planes[0];

	magI += Scalar::all(1);
	log(magI, magI);                //ת���������߶�(logarithmic scale)

									//����������л��У����Ƶ�׽��вü�
	magI = magI(Rect(0, 0, magI.cols&-2, magI.rows&-2));

	//�������и���Ҷͼ���е����ޣ�ʹ��ԭ��λ��ͼ������
	int cx = magI.cols / 2;
	int cy = magI.rows / 2;

	Mat q0(magI, Rect(0, 0, cx, cy));       //���Ͻ�ͼ�񻮶�ROI����
	Mat q1(magI, Rect(cx, 0, cx, cy));      //���Ͻ�ͼ��
	Mat q2(magI, Rect(0, cy, cx, cy));      //���½�ͼ��
	Mat q3(magI, Rect(cx, cy, cx, cy));     //���½�ͼ��

											//�任���ϽǺ����½�����
	Mat tmp;
	q0.copyTo(tmp);
	q3.copyTo(q0);
	tmp.copyTo(q3);

	//�任���ϽǺ����½�����
	q1.copyTo(tmp);
	q2.copyTo(q1);
	tmp.copyTo(q2);

	//��һ��������0-1֮��ĸ�����������任Ϊ���ӵ�ͼ���ʽ
	normalize(magI, magI, 0, 1, CV_MINMAX);

	outfile1 << magI << endl;
	outfile1.close();
	imshow("����ͼ��", I);
	imshow("Ƶ��ͼ", magI);
	waitKey(0);


	return 0;
}