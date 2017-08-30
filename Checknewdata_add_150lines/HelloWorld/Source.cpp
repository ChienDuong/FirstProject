// da chuyen sang generate data_ Vertical  : truyen lai vao bien generate data khi xuat ra
// chi can hoan thanh initial model nua cho vertical
#include <iostream>
#include <fstream>
#include <vector> // giup tao dc array ma co dimensional cuc cao
#include <time.h>
#include <stdlib.h>// dung de goi ham malloc va calloc trong C
#include<iomanip>
#include<cmath> 
#include <windows.h>

//#include<iostream>
#include<opencv2\highgui\highgui.hpp>
#include<opencv2\imgproc.hpp>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"

#include "opencv2/calib3d/calib3d.hpp"
//#include "opencv2/nonfree/nonfree.hpp"
using namespace cv;
using namespace std;
//using std::vector;
const double spacing_min = 0.04;
const double spacing_max = 0.073;
const double ConstraintsX = 0.065;
const double dist_LL_Model = 0.03;
int switch_V_H = 2;// 0 : vertical, 1 : horizol, 2: both // hor o tren, ver o duoi

void StartProgram1(string filename, string filename1)
{
	ShellExecuteA(NULL, "open" , filename.c_str(), filename1.c_str(), NULL, SW_SHOW);

}

void GetPointsForWarp(vector<cv::Point2f> &pts)
{
	pts.push_back(cv::Point2f(0.0, 0.0));
	pts.push_back(cv::Point2f(200.041, 0.0581395));
	pts.push_back(cv::Point2f(-0.0639535, 200.006));
	pts.push_back(cv::Point2f(199.977, 200.07));

}
void GetPointsCADForWarp(vector<cv::Point2f> &pts)
{
	pts.push_back(cv::Point2f(219.3799, 151.2251));
	pts.push_back(cv::Point2f(419.3799, 151.2251));
	pts.push_back(cv::Point2f(219.3799, 351.2251));
	pts.push_back(cv::Point2f(419.3799, 351.2251));

}


double dist_P2P(vector<double> &X1, vector<double> &X2)
{
	double d = sqrt(pow((X1[0] - X2[0]), 2) + pow((X1[1] - X2[1]), 2));
	return d;
};
double dist_P2L(vector<double> &X1, vector<double> &X2)
{
	double dl = (X2.at(0) * X1.at(0) - X1.at(1) + X2.at(1)) / sqrt(pow(X2.at(0), 2) + pow(X2.at(1), 2));
	return dl;
};

// input, toa do diem, va bac cua phuong trinh
vector <double> Polyfit(vector<double> &X_point, vector<double> &Y_point, int n)
{
	int i, j, k, N;
	N = X_point.size();
	cout.precision(4);                        //set precision
	cout.setf(ios::fixed);
	vector<double> X(2 * n + 1);                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	for (i = 0; i < 2 * n + 1; i++)
	{
		X[i] = 0;
		for (j = 0; j < N; j++)
			X[i] = X[i] + pow(X_point.at(j), i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
	}
	vector<double> a(n + 1);
	vector<double> am(n + 1);// them vao de cho dung thu tu
	vector<vector<double>>B(n + 1, vector<double>(n + 2));
	//double B[n + 1][n + 2], a[n + 1];            //B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients
	for (i = 0; i <= n; i++)
	for (j = 0; j <= n; j++)
		B[i][j] = X[i + j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix
	vector<double> Y(n + 1);
	//double Y[n + 1];                    //Array to store the values of sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	for (i = 0; i < n + 1; i++)
	{
		Y[i] = 0;
		for (j = 0; j < N; j++)
			Y[i] = Y[i] + pow(X_point.at(j), i)*Y_point.at(j);        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	}
	for (i = 0; i <= n; i++)
		B[i][n + 1] = Y[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)
	n = n + 1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations

	for (i = 0; i < n; i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
	for (k = i + 1; k < n; k++)
	if (B[i][i] < B[k][i])
	for (j = 0; j <= n; j++)
	{
		double temp = B[i][j];
		B[i][j] = B[k][j];
		B[k][j] = temp;
	}

	for (i = 0; i < n - 1; i++)            //loop to perform the gauss elimination
	for (k = i + 1; k < n; k++)
	{
		double t = B[k][i] / B[i][i];
		for (j = 0; j <= n; j++)
			B[k][j] = B[k][j] - t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
	}
	for (i = n - 1; i >= 0; i--)                //back-substitution
	{                        //x is an array whose values correspond to the values of x,y,z..
		a[i] = B[i][n];                //make the variable to be calculated equal to the rhs of the last equation
		for (j = 0; j < n; j++)
		if (j != i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
			a[i] = a[i] - B[i][j] * a[j];
		a[i] = a[i] / B[i][i];            //now finally divide the rhs by the coefficient of the variable to be calculated
	}

	am = a;
	for (int i = 0; i < n; i++)
	{
		a[i] = am[n - i - 1];
	}

	return a;
};
void Ref_model(vector<double>  &s, vector<int> &pos_outlier, vector<vector<int>> & dataXY, vector<vector<double>> &generate_data, vector<vector<double>> &XY, double ConstraintsX)
{

	vector <double> predict_point(2, 0);
	int range_of_candidate = 3;
	vector<double> current_point(2);
	current_point = XY[dataXY[pos_outlier.at(1)][pos_outlier.at(0)]]; \

		// Check If this is a first point, so just predict by +0.065  to X_axes do not need ref_model
	if (pos_outlier.at(1) == 1)
	{
		predict_point = current_point;
		predict_point.at(0) == current_point.at(0) + ConstraintsX;
		exit;
	}
	vector<double> ref_candidate(range_of_candidate * 2 + 1);// luu cac gia tri ref theo x
	if (pos_outlier.at(0) <= (generate_data[0].size() - range_of_candidate))
	{
		for (int i = -range_of_candidate; i < range_of_candidate + 1; i++)
		{
			ref_candidate.at(i + range_of_candidate) = XY[generate_data[pos_outlier.at(1) - 1][pos_outlier.at(0) + i]][0];// OK chay tu 0 den 6, mat lab chay tu 1 den 7
		}
	}
	else
	{
		ref_candidate[0] = XY[generate_data[pos_outlier.at(1) - 1][pos_outlier.at(0)]][0];
	}
	double predict_value;
	double min_predict, min_predict_temp;
	int predict_ref;// cai nay se nho hon trong matlab 1
	min_predict = abs(ref_candidate.at(0) - current_point.at(0));
	// Choose the below 7 points, and choose the one have smallest different X
	for (int i = 0; i < range_of_candidate * 2 + 1; i++)
	{
		//ref_candidate.at(0);
		min_predict_temp = abs(ref_candidate.at(i) - current_point.at(0));
		if (min_predict>min_predict_temp)
		{
			min_predict = min_predict_temp;
			predict_ref = i;
		}

	}
	int Ref_indexX, Ref_indexY;
	Ref_indexX = pos_outlier.at(0) - range_of_candidate + predict_ref;// do pos_outlier giam 1, range_of_candidate giam 1 => o day ko tru 1 nua so voi matlab
	Ref_indexY = pos_outlier.at(1) - 1;
	vector<double> Ref_point(2);
	vector<double>Ref_nextpoint(2);
	vector<double>predict_vector(2);
	Ref_point = XY[generate_data[Ref_indexY][Ref_indexX]]; // XY[300] tra ve toa do diem 300 do luon
	Ref_nextpoint = XY[generate_data[Ref_indexY][Ref_indexX + 1]];
	predict_vector = { Ref_nextpoint.at(0) - Ref_point.at(0), Ref_nextpoint.at(1) - Ref_point.at(1) };
	// Ref_vector = diffrent with next and current reference points
	if ((predict_vector.at(0) + predict_vector.at(1)) == 0)
	{
		Ref_nextpoint = XY[generate_data[Ref_indexY][Ref_indexX + 2]];
		predict_vector = { Ref_nextpoint.at(0) - Ref_point.at(0), Ref_nextpoint.at(1) - Ref_point.at(1) };
	}
	vector<double> predict_point2(2);
	vector<double> predict_point1(2);
	predict_point2 = current_point;
	predict_point2.at(0) = current_point.at(0) + ConstraintsX;
	predict_point1 = { current_point.at(0) + predict_vector.at(0), current_point.at(1) + predict_vector.at(1) };
	predict_point = { (predict_point1.at(0) + predict_point2.at(0)) / 2, (predict_point1.at(1) + predict_point2.at(1)) / 2 };
	s.at(0) = (predict_point.at(0));
	s.at(1) = (predict_point.at(1));
}
void InitialRef_model(vector<double>  &s, vector<int> &pos_outlier, vector<vector<int>> & dataXY, vector<vector<double>> &generate_data, vector<vector<double>> &XY, double ConstraintsX)
{

	vector <double> predict_point(2, 0);
	int range_of_candidate = 3;
	vector<double> current_point(2);
	current_point = XY[dataXY[pos_outlier.at(1)][pos_outlier.at(0)]]; \

	
	vector<double> ref_candidate(range_of_candidate * 2 + 1);// luu cac gia tri ref theo x
	if (pos_outlier.at(0) <= (generate_data[0].size() - range_of_candidate) && pos_outlier.at(0)>range_of_candidate)
	{
	//	vector<double> ref_candidate(range_of_candidate * 2 + 1);// luu cac gia tri ref theo x
		for (int i = -range_of_candidate; i < range_of_candidate + 1; i++)
		{
			ref_candidate.at(i + range_of_candidate) = XY[generate_data[pos_outlier.at(1) + 1][pos_outlier.at(0) + i]][0];// OK chay tu 0 den 6, mat lab chay tu 1 den 7
		}
	}
	else
	{
		range_of_candidate = 0;
		ref_candidate.resize(range_of_candidate * 2 + 1);// luu cac gia tri ref theo x
		ref_candidate[0] = XY[generate_data[pos_outlier.at(1) + 1][pos_outlier.at(0)]][0];
		
	}
	double predict_value;
	double min_predict, min_predict_temp;
	int predict_ref;// cai nay se nho hon trong matlab 1
	min_predict = abs(ref_candidate.at(0) - current_point.at(0));
	// Choose the below 7 points, and choose the one have smallest different X
	for (int i = 0; i < range_of_candidate * 2 + 1; i++)
	{
		//ref_candidate.at(0);
		min_predict_temp = abs(ref_candidate.at(i) - current_point.at(0));
		if (min_predict>min_predict_temp)
		{
			min_predict = min_predict_temp;
			predict_ref = i;
		}

	}
	int Ref_indexX, Ref_indexY;
	Ref_indexX = pos_outlier.at(0) - range_of_candidate + predict_ref;// do pos_outlier giam 1, range_of_candidate giam 1 => o day ko tru 1 nua so voi matlab
	Ref_indexY = pos_outlier.at(1) + 1;
	vector<double> Ref_point(2);
	vector<double>Ref_nextpoint(2);
	vector<double>predict_vector(2);
	Ref_point = XY[generate_data[Ref_indexY][Ref_indexX]]; // XY[300] tra ve toa do diem 300 do luon
	Ref_nextpoint = XY[generate_data[Ref_indexY][Ref_indexX + 1]];
	predict_vector = { Ref_nextpoint.at(0) - Ref_point.at(0), Ref_nextpoint.at(1) - Ref_point.at(1) };
	// Ref_vector = diffrent with next and current reference points
	if ((predict_vector.at(0) + predict_vector.at(1)) == 0)
	{
		Ref_nextpoint = XY[generate_data[Ref_indexY][Ref_indexX + 2]];
		predict_vector = { Ref_nextpoint.at(0) - Ref_point.at(0), Ref_nextpoint.at(1) - Ref_point.at(1) };
	}
	vector<double> predict_point2(2);
	vector<double> predict_point1(2);
	predict_point2 = current_point;
	predict_point2.at(0) = current_point.at(0) + ConstraintsX;
	predict_point1 = { current_point.at(0) + predict_vector.at(0), current_point.at(1) + predict_vector.at(1) };
	predict_point = { (predict_point1.at(0) + predict_point2.at(0)) / 2, (predict_point1.at(1) + predict_point2.at(1)) / 2 };
	s.at(0) = (predict_point.at(0));
	s.at(1) = (predict_point.at(1));
}

//  return the number of column and rows
int * readFile(char* filename, int* a, int*b)
{

	FILE* fp;
	int countLines = 0;
	int i;
	int columns = 0;
	int maxcolumns = 0;
	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("We can't open the file.");
		fclose(fp);
		exit(1);
	}
	else
	{
		int n;
		std::ifstream inFile(filename);
		n = std::count(std::istreambuf_iterator<char>(inFile),
			std::istreambuf_iterator<char>(), '\n');
		//static int r[200];
		static int *r = (int *)malloc(n * sizeof(int));
		//static std::vector<int> r(n);

		while ((i = fgetc(fp)) != EOF)
		{

			// lay cai co hang lon nhat lam chuan
			if (i == ',' || i == ' ')
			{
				++columns;
			}

			if (i == '\n')
			{
				r[countLines] = columns;
				if (maxcolumns < columns)
				{
					maxcolumns = columns;
					columns = 0;
				}
				else
				{
					columns = 0;  // reset for next line.
				}

				countLines++;
			}

		}

		printf(">Numrows filename.txt: %d\n", countLines);
		printf(">Num columns a.txt: %d\n", (maxcolumns + 1));
		*a = (countLines);
		*b = maxcolumns + 1;
		return r;
	}
}

int main()
{

	
	clock_t tic, toc;
	tic = clock();
	cout << "tic:" << tic << endl;
	vector<cv::Point2f> SCAN (4);
	ifstream PointsForWarp("GetPointsForWarp.TXT");//105L.H.G1X.TXT

	if (PointsForWarp.is_open())
	{
		// cach muon tao 1 mang
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (j==0)
				{
					PointsForWarp >> SCAN[i].x;
				}
				else
				{
					PointsForWarp >> SCAN[i].y;
				}
								
			}


		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}
	PointsForWarp.close();
	
	//GetPointsForWarp(SCAN);
	/*vector<cv::Point2f> CAD;
	GetPointsCADForWarp(CAD);*/
	vector<cv::Point2f> CAD(4);
	ifstream PointsCADForWarp("GetPointsCADForWarp.TXT");//105L.H.G1X.TXT

	if (PointsCADForWarp.is_open())
	{
		// cach muon tao 1 mang
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				if (j == 0)
				{
					PointsCADForWarp >> CAD[i].x;
				}
				else
				{
					PointsCADForWarp >> CAD[i].y;
				}

			}


		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}
	PointsCADForWarp.close();
	Mat H = findHomography(SCAN, CAD);

	cout << H;


	int a, b;
	int a_H, b_H;
	int *p;
	//**************************Read data******************************************************************************************************//
	toc = clock();
	printf("time to homography: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	// Horizoltal data
	p = readFile("testdata105L.H.G2X_150.TXT", &a_H, &b_H);//105L.H.G1X.TXT
	vector<vector<double> > Hor_X(a_H, vector<double>(b_H));
	toc = clock();
	printf("read number of rows,columns: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	// Đọc các giá trị của ma trận từ file
	ifstream H_H_L("testdata105L.H.G2X_150.TXT");//105L.H.G1X.TXT

	if (H_H_L.is_open())
	{
		// cach muon tao 1 mang
		for (int i = 0; i < a_H; i++)
		{
			for (int j = 0; j < b_H; j++)
			{
				// xoa nhung gia tri lon hon so column cua dong do
				if (j> *(p + i))
				{

					Hor_X[i][j] = 0;
					continue;
				}
				H_H_L >> Hor_X[i][j];				
			}


		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}

	cout << "Write succesfully" << endl;
	
	H_H_L.close();
	toc = clock();
	printf("Store Horizoltal_X data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		
	// Luu gia tri Y*****************************************************************************//
	vector<vector<double> > Hor_Y(a_H, vector<double>(b_H));
	// Đọc các giá trị của ma trận từ file
	ifstream H_V_L("testdata105L.H.G2Y_150.TXT");//105L.H.G1Y.TXT

	if (H_V_L.is_open())
	{
		
		for (int i = 0; i < a_H; i++)// a=105 b= 3138
		{
			for (int j = 0; j < b_H; j++)
			{
				// xoa nhung gia tri lon hon so column cua dong do
				if (j> *(p + i))
				{

					Hor_Y[i][j] = 0;
					continue;
				}
				H_V_L >> Hor_Y[i][j];			
			}


		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}
	H_V_L.close();
	toc = clock();
	printf("Store Horizoltal_Y data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		
	int t_H = a_H*b_H;

	vector<vector<double> > Hor_XY(t_H, vector<double>(2));

	for (int i = 0; i < b_H; i++)//a=105 b=3138, gia tri cua Ver_XY dc in theo hang cua Ver_Y va Ver_X
	{
		for (int j = 0; j < a_H; j++)//105
		{
			Hor_XY[j + a_H*i][0] = Hor_X[j][i];  //0 la X thi chinh o day
			Hor_XY[j + a_H*i][1] = Hor_Y[j][i];
		}
	}
	toc = clock();
	printf("Store Horizoltal data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	//tic = clock();
	// Vertical data
	p = readFile("testdata105L.V.G2X_150.txt", &a, &b);//105L.V.G1X.TXT
	vector<vector<double> > Ver_X(a, vector<double>(b));
	toc = clock();
	printf("read number of rows,columns of Vertical: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	// Đọc các giá trị của ma trận từ file
	ifstream H_L("testdata105L.V.G2X_150.TXT");//105L.V.G1X.TXT

	if (H_L.is_open())
	{
		// cach muon tao 1 mang
		for (int i = 0; i < a; i++)
		{
			for (int j = 0; j < b; j++)
			{
				// xoa nhung gia tri lon hon so column cua dong do
				if (j> *(p + i))
				{

					Ver_X[i][j] = 0;
					continue;
				}
				H_L >> Ver_X[i][j];// su dung ham nay phai chuyen dau phay sang khoang cach				
			}

		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}
	H_L.close();
	toc = clock();
	printf("Store Vertical_X data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	// Luu gia tri Y*****************************************************************************//
	vector<vector<double> > Ver_Y(a, vector<double>(b));
	// Đọc các giá trị của ma trận từ file
	ifstream V_L("testdata105L.V.G2Y_150.TXT");//105L.V.G1Y.TXT

	if (V_L.is_open())
	{
		// cach muon tao 1 mang
		for (int i = 0; i < a; i++)// a=105 b= 3138
		{
			for (int j = 0; j < b; j++)
			{
				// xoa nhung gia tri lon hon so column cua dong do
				if (j> *(p + i))
				{

					Ver_Y[i][j] = 0;
					continue;
				}
				V_L >> Ver_Y[i][j];// su dung ham nay phai chuyen dau phay sang khoang cach				
			}


		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}
	V_L.close();
	toc = clock();
	printf("Store Vertical_Y data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	int t = a*b;

	vector<vector<double> > Ver_XY(t, vector<double>(2));

	for (int i = 0; i < b; i++)//a=105 b=3138, gia tri cua Ver_XY dc in theo hang cua Ver_Y va Ver_X
	{
		for (int j = 0; j < a; j++)//105
		{
			Ver_XY[j + a*i][0] = Ver_X[j][i];
			Ver_XY[j + a*i][1] = Ver_Y[j][i];
		}
	}
	toc = clock();
	printf("Stor Vertical data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);


	//1************************** finished Read data******************************************************************************************************//

	//2***********************************Creat New table******************************************************************************************************//
	// Horizoltal data
	vector<vector<int>>Hor_dataXY(a_H, vector < int >(b_H)); // 105x3138
	int index_H = 0;
	for (int j = 0; j < Hor_X[0].size(); j++)// chinh la b_H
	{
		for (int i = 0; i < Hor_X.size(); i++)// chinh la a_H
		{
			Hor_dataXY[i][j] = index_H;
			index_H = index_H + 1;

		}
	}

	toc = clock();
	printf("Creat Newtable H in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	// Vertical data
	// NOTE: Created tranposed matrix
	vector<vector<int>>Ver_dataXY(a, vector < int >(b)); // 105x3138
	int index = 0;
	for (int j = 0; j < Ver_X[0].size(); j++)// chinh la b
	{
		for (int i = 0; i < Ver_X.size(); i++)// chinh la a
		{
			Ver_dataXY[i][j] = index;
			index = index + 1;

		}
	}
	toc = clock();
	printf("Creat new table Y in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

	//2*********************************** Finihsed New table******************************************************************************************************//
	
	//4***********************************Plot Horizontal boundary up and down *********************************************************************************//
	//tic = clock();
	cout << "Plot boundary up and down" << endl;
	//vector<vector<double>>cross_point_UandD(1, vector<double>(2)) ;
	vector<vector<double>>cross_point_UandD;
	for (int indexX_2 = 0; indexX_2 < Hor_dataXY.size(); indexX_2 = (indexX_2 + Hor_dataXY.size() - 1))   // !!X_2 HOr = 105x3138 , X_1 Ver =105x3111
	{
		int firstY_1 = 0;// store the index of Vertical points of first cross point
		int firstpreviousY_2 = 0;
		int secondpreviousY_2 = 0;
		int tempY_2 = 0;// store the average number of Horizoltal points between 2 vertical lines 
		int indexX_1 = 0;
		for (int indexY_2 = 0; indexY_2 < Hor_dataXY[0].size(); indexY_2++)// cua horizoltal 3138,6317 xet tat ca cac diem cua current horizoltal
		{
			for (int indexY_1 = 0; indexY_1 < Ver_dataXY[0].size(); indexY_1++)// cua vertical 3111, 6223
			{
				if (indexX_1>Ver_dataXY.size() - 1)//!! Chu y do khac nhau giua matLab va C nen phai tru 1, ngoai line 105 bo?
				{
					break;
				}
				if (indexX_1==0)
				{
					double cross = dist_P2P(Hor_XY[Hor_dataXY[indexX_2][indexY_2]], Ver_XY[Ver_dataXY[indexX_1][indexY_1]]); // xet tung diem cua horizoltal voi vertical
					if (cross == 0)
					{
						cross_point_UandD.push_back(vector<double>{Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(0), Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(1)});// ko the push back 1 vector 2 mang vao 2 mang, ma hay nhet diem trong mang 2 chieu do (tuong ung vs 1 vector)
						firstY_1 = indexY_1;
						firstpreviousY_2 = indexY_2;						
						indexX_1 = indexX_1 + 1;// xet hang moi cua vertical
						break;// xet hang moi cua vertical
					}
					else
					{
						continue;// xet diem moi trong hang cua vertical
					}
					
				}
				if (indexX_1 == 1)
				{
					double cross = dist_P2P(Hor_XY[Hor_dataXY[indexX_2][indexY_2]], Ver_XY[Ver_dataXY[indexX_1][indexY_1]]); // xet tung diem cua horizoltal voi vertical
					if (cross == 0)
					{
						cross_point_UandD.push_back(vector<double>{Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(0), Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(1)});// ko the push back 1 vector 2 mang vao 2 mang, ma hay nhet diem trong mang 2 chieu do (tuong ung vs 1 vector)
						firstY_1 = indexY_1;
						secondpreviousY_2 = indexY_2;
						indexX_1 = indexX_1 + 1;
						break;
					}
					else
					{
						continue;// xet diem moi trong hang cua vertical
					}

				}
				tempY_2 = secondpreviousY_2 - firstpreviousY_2;
			
				
				//double cross = dist_P2P(Hor_XY[Hor_dataXY[indexX_2][indexY_2]], Ver_XY[Ver_dataXY[indexX_1][indexY_1]]); // xet tung diem cua horizoltal voi vertical
				//if (cross == 0)
				//{
				//	cross_point_UandD.push_back(vector<double>{Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(0), Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(1)});// ko the push back 1 vector 2 mang vao 2 mang, ma hay nhet diem trong mang 2 chieu do (tuong ung vs 1 vector)
				//	indexX_1 = indexX_1 + 1;
				//	indexY_2 = indexY_2 + tempY_2-2;
				//	break;
				//}
			}
			
		}
		for (int indexY_2 = 0; indexY_2 < Hor_dataXY[0].size(); indexY_2++)// cua horizoltal 3138, xet tat ca cac diem cua current horizoltal, Y_2: horizoltal, Y_1, vertical
		{
			for (int indexY_1 = firstY_1 - 2; indexY_1 < firstY_1 + 2; indexY_1++)// bay h chi xet 5 diem dao dong xung quanh, diem dau tien cross
			{
				if (indexX_1>Ver_dataXY.size() - 1)//!! Chu y do khac nhau giua matLab va C nen phai tru 1, ngoai line 105 bo?
				{
					break;
				}
						
				if (firstY_1<0) // check the cross is the first point or not!
				{
					cout << "Warning!!!!!!!!!!!!!!!!!!!!, can not find cross point. Check input data again";
					continue;
				}
				double cross = dist_P2P(Hor_XY[Hor_dataXY[indexX_2][indexY_2]], Ver_XY[Ver_dataXY[indexX_1][indexY_1]]); // xet tung diem cua horizoltal voi vertical
				if (cross == 0)
				{
					cross_point_UandD.push_back(vector<double>{Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(0), Hor_XY[Hor_dataXY[indexX_2][indexY_2]].at(1)});// ko the push back 1 vector 2 mang vao 2 mang, ma hay nhet diem trong mang 2 chieu do (tuong ung vs 1 vector)
					indexX_1 = indexX_1 + 1;
					indexY_2 = indexY_2 + tempY_2 - 2;
					break;
				}
			}
		}
	}
	toc = clock();
	printf("Boundary up and down in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	////4***********************************Finished Plot Horizontal boundary up and down ********************************************************************************//
	//4***********************************Plot Vertical boundary left and right *********************************************************************************//
	cout << "Plot Vertical boundary left and right" << endl;
	//vector<vector<double>>cross_point_UandD(1, vector<double>(2)) ;
	vector<vector<double>>cross_point_LandR;
	for (int indexX_2 = 0; indexX_2 < Ver_dataXY.size(); indexX_2 = (indexX_2 + Ver_dataXY.size() - 1))   // !! HOr = 105x3138 , Ver =105x3111
	{
		int firstY_1 = 0;// store the index of Vertical points of first cross point
		int firstpreviousY_2 = 0;
		int secondpreviousY_2 = 0;
		int tempY_2 = 0;// store the average number of Horizoltal points between 2 vertical lines 
		int indexX_1 = 0;
		for (int indexY_2 = 0; indexY_2 < Ver_dataXY[0].size(); indexY_2++)// cua horizoltal 3138, xet tat ca cac diem cua current horizoltal
		{
			for (int indexY_1 = 0; indexY_1 < Hor_dataXY[0].size(); indexY_1++)// cua vertical 3111
			{
				if (indexX_1>Hor_dataXY.size() - 1)//!! Chu y do khac nhau giua matLab va C nen phai tru 1, ngoai line 105 bo?
				{
					break;
				}
				if (indexX_1 == 0)
				{
					double cross = dist_P2P(Ver_XY[Ver_dataXY[indexX_2][indexY_2]], Hor_XY[Hor_dataXY[indexX_1][indexY_1]]); // xet tung diem cua horizoltal voi vertical
					if (cross == 0)
					{
						cross_point_LandR.push_back(vector<double>{Ver_XY[Ver_dataXY[indexX_2][indexY_2]].at(0), Ver_XY[Ver_dataXY[indexX_2][indexY_2]].at(1)});// ko the push back 1 vector 2 mang vao 2 mang, ma hay nhet diem trong mang 2 chieu do (tuong ung vs 1 vector)
						firstY_1 = indexY_1;
						firstpreviousY_2 = indexY_2;
						indexX_1 = indexX_1 + 1;
						break;// xet hang moi cua vertical
					}
					else
					{
						continue;// xet diem moi trong hang cua vertical
					}

				}
				if (indexX_1 == 1)
				{
					double cross = dist_P2P(Ver_XY[Ver_dataXY[indexX_2][indexY_2]], Hor_XY[Hor_dataXY[indexX_1][indexY_1]]); // xet tung diem cua horizoltal voi vertical
					if (cross == 0)
					{
						cross_point_LandR.push_back(vector<double>{Ver_XY[Ver_dataXY[indexX_2][indexY_2]].at(0), Ver_XY[Ver_dataXY[indexX_2][indexY_2]].at(1)});// ko the push back 1 vector 2 mang vao 2 mang, ma hay nhet diem trong mang 2 chieu do (tuong ung vs 1 vector)
						firstY_1 = indexY_1;
						secondpreviousY_2 = indexY_2;
						indexX_1 = indexX_1 + 1;
						break;
					}
					else
					{
						continue;// xet diem moi trong hang cua vertical
					}

				}
				tempY_2 = secondpreviousY_2 - firstpreviousY_2;

			}

		}
		for (int indexY_2 = 0; indexY_2 < Ver_dataXY[0].size(); indexY_2++)// cua horizoltal 3138, xet tat ca cac diem cua current horizoltal
		{
			for (int indexY_1 = firstY_1 - 2; indexY_1 < firstY_1 + 2; indexY_1++)// bay h chi xet 5 diem dao dong xung quanh, diem dau tien cross
			{
				if (indexX_1>Hor_dataXY.size() - 1)//!! Chu y do khac nhau giua matLab va C nen phai tru 1, ngoai line 105 bo?
				{
					break;
				}

				
				double cross = dist_P2P(Ver_XY[Ver_dataXY[indexX_2][indexY_2]], Hor_XY[Hor_dataXY[indexX_1][indexY_1]]); // xet tung diem cua 1:horizoltal voi 2:vertical
				if (cross == 0)
				{
					cross_point_LandR.push_back(vector<double>{Ver_XY[Ver_dataXY[indexX_2][indexY_2]].at(0), Ver_XY[Ver_dataXY[indexX_2][indexY_2]].at(1)});// ko the push back 1 vector 2 mang vao 2 mang, ma hay nhet diem trong mang 2 chieu do (tuong ung vs 1 vector)
					indexX_1 = indexX_1 + 1;
					indexY_2 = indexY_2 + tempY_2 - 2;
					break;
				}
			}
		}
	}
	toc = clock();
	printf("Boundary left and right: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	/*toc = clock();
	printf("Boundary in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);*/
	
//4***********************************Finished Plot Vertical boundary left and right ********************************************************************************//
	//5***********************************Process vertical select line ********************************************************************************//
	
	if (switch_V_H == 1 || switch_V_H == 2)
	{

	

		//5***********************************Process Hor select line ********************************************************************************//
		//6*****************************************Initial first line*******************************************************************************************//
		/*tic = clock();*/
		vector<vector<int>> table_outlier; // khoi tao = 1x1
		vector<int> outlier_flag;// = { 0 };
		int indexX;
		int indexY;
		int firstindex;
		int initial_model = 0;
		int initial_flag = 0;
		int firstY = -1;
		int have_predict = 0;
		int i;
		int index1;
		int index2;
		vector<double> predict_point(2);
		predict_point.reserve(2);
		vector<int>sliding_window(2);
		vector<int> pos_outlier(2);
		vector <double> LL_model = { 1, 2 };// linear line model tim bang polyfit sua sau
		vector<vector<double> > generate_Hor_data(Hor_dataXY.size(), vector<double>(Hor_dataXY[0].size()));
		while (1)
		{
			//tic = clock();
			if (initial_flag == 1)
			{
				firstY = firstY - 1;
			}
			else
			{
				firstY = firstY + 1;
			}
			if (firstY == -1)
			{
				break;
			}
			cout << "firstY =" << firstY << endl;
			indexX = 0;
			indexY = firstY;
			//outlier_flag[firstY]=0; // cai nay nen dung pushback
			outlier_flag.push_back(0);
			if (Hor_dataXY.size() == 0)
			{
				continue;
			}


			firstindex = Hor_dataXY[indexY][indexX];
			sliding_window = { (indexX + 1), indexY };
			vector<int> Line;
			Line.push_back(firstindex);  // Line ban dau =1, diem dau tien trong Hor_dataXY, sau do move sang diem moi

			i = 0;
		
			//	// check end line!
			while (1)
			{
				if (outlier_flag.back() == 1 && initial_flag == 0)
				{
					break;
				}
				i = i + 1;
				if (i == Hor_dataXY[0].size())
				{
					break;
				}

				//		// Make a Line from 5 points by optimization projection errors
				if (!Line.empty())
				{
					index1 = Line.back();// -1so voi matlab, Vi line luu giu size cua Hor_XY, nhung ma chi so xuat phat tu 0 => -1
				}
				sliding_window = { indexX + i, indexY };// o day 6 0 thi matlab la 7 1
				index2 = Hor_dataXY[sliding_window[1]][sliding_window[0]];// nho hon matlab 1 so
				
				if (have_predict == 1)
				{
					if (dist_P2P(predict_point, Hor_XY[index2]) < spacing_min)
					{
						Hor_XY[index2] = predict_point;
						if (Line.back() != index2)
						{
							Line.push_back(index2);
						}
						have_predict = 0;
						continue;

					}
					else
					{
						if (Hor_XY[index2][0] < predict_point.at(1))
						{
							Hor_XY[index2] = predict_point;
							if (Line.back() != index2)
							{
								Line.push_back(index2);
							}
							have_predict = 0;
							continue;
						}
						//clock_t tic, toc;
						//	tic = clock();
						Hor_XY.push_back(predict_point);// truong hop mat diem, thi se tao them 1 diem luu vao cuoi cua Hor_XY, va index cua no dc luu vao line, tu do de ve thi ve cac index chua trong line
						Line.push_back(Hor_XY[0].size());
						have_predict = 0;
						//			toc = clock();
						//		printf("Loading file aa: %f seconds", (double)(toc - tic) / CLOCKS_PER_SEC);
						continue;
					}

				}
				double d = dist_P2P(Hor_XY[index1], Hor_XY[index2]);// distance chinh xac
				//		//TESTCASE
				//		//cout << d; // 0.064
				if (d < spacing_min)
				{
					cout << "Warning : Constraints1(CT1): 0.04<distance" << endl;
					Hor_XY[index2] = Hor_XY[index1];
					if (((Hor_XY[index2][0]) + (Hor_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					continue;
				}
				else if (d > spacing_max) {
					cout << "Warning : Constraints1(CT1): distance<0.073" << endl;
					/*cout << Hor_XY[index1 - 105].at(0) << " ";
					cout << Hor_XY[index1 - 105].at(1) << endl;
					cout << Hor_XY[index1].at(0)<<" ";
					cout << Hor_XY[index1].at(1)<<endl;
					
					cout << Hor_XY[index2].at(0)<<" ";
					cout << Hor_XY[index2].at(1)*/;
					//cout << Hor_XY[index2];*/
					if (((Hor_XY[index2][0]) + (Hor_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					else
					{
						HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
						SetConsoleTextAttribute(h, FOREGROUND_RED | FOREGROUND_INTENSITY);
						cout << "Warning: outlier" << endl;
						//                     % generate the points and run the loop again
						//                  % can avoid two kinds of situation the outliers point
						SetConsoleTextAttribute(h, 15);
						//	outlier_flag(firstY) = 1;
						outlier_flag.back() = 1;
						//---------------------------------------------predict point----------------------------------------------------//
						if (initial_flag == 1)
						{


							pos_outlier = { i - 1, firstY };// phai tru 1 do xuat phat tu 0, ma i o day khoi tao la 1
							table_outlier.push_back(pos_outlier);
							
							InitialRef_model(predict_point, pos_outlier, Hor_dataXY, generate_Hor_data, Hor_XY, ConstraintsX);
							
							have_predict = 1;
							//toc = clock();
							//					printf("Ref model outside: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

							i = i - 1;
						}
						//---------------------------------------------predict point----------------------------------------------------//
					}
					continue;

				}
				else if (Hor_XY[index2][0] < Hor_XY[index1][0])
				{
					cout << Hor_XY[index2][0];
					cout << Hor_XY[index1][0];
					cout << "Warning : Next point should in right side" << endl;
					continue;
				}

				if (Line.size() == 5)
				{
					vector<double>point_X(5);
					vector<double>point_Y(5);
					for (int i = 0; i < 5; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X[i] = Hor_XY[Line.at(i)][0];// test case: -0.157,-0.093,-0.0291,0.0349,0.0988
						point_Y[i] = Hor_XY[Line.at(i)][1];//test case: 5 so 23.1686

					}
					LL_model = Polyfit(point_X, point_Y, 1); // test case: 0 23.1686
				}
				else if (Line.size() > 5 && Line.size() <= 10)
				{
					if (dist_P2L(Hor_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						//outlier_flag(firstY) = 1;
						outlier_flag.back() = 1;
						if (initial_flag == 1)
						{


							// Predict points//
							pos_outlier = { i - 1, firstY };
							table_outlier.push_back(pos_outlier);
							InitialRef_model(predict_point, pos_outlier, Hor_dataXY, generate_Hor_data, Hor_XY, ConstraintsX);
							have_predict = 1;
							i = i - 1;
						}
						continue;
					};
					vector<double>point_X1;
					vector<double>point_Y1;  // moi vong lap deu ve lai

					for (int i = 0; i < Line.size(); i++) // Line tu 0 den 4, mat lab la 1 den 5
					{
						//Hor_XY[300];
						point_X1.push_back(Hor_XY[Line.at(i)][0]);
						point_Y1.push_back(Hor_XY[Line.at(i)][1]);
					}

					LL_model = Polyfit(point_X1, point_Y1, 1);

				}

				else if (Line.size() > 10)
				{
					if (dist_P2L(Hor_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						//	outlier_flag(firstY) = 1;
						outlier_flag.back() = 1;
						if (initial_flag == 1)
						{


							// Predict points//
							pos_outlier = { i - 1, firstY };
							table_outlier.push_back(pos_outlier);
							InitialRef_model(predict_point, pos_outlier, Hor_dataXY, generate_Hor_data, Hor_XY, ConstraintsX);
							have_predict = 1;
							i = i - 1;
						}
						continue;
					}
					vector<double>point_X2(10);
					vector<double>point_Y2(10);
					for (int i = (Line.size() - 9); i < Line.size() + 1; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X2.at(i - (Line.size() - 9)) = Hor_XY[Line.at(i - 1)][0];
						point_Y2.at(i - (Line.size() - 9)) = Hor_XY[Line.at(i - 1)][1];
					}
					LL_model = Polyfit(point_X2, point_Y2, 1);
				}
				if (Line.back() != index2)
				{
					Line.push_back(index2);
				}

			}// vong lap while thu nhat chay tat ca cac diem tren 1 hang
			for (int i = 0; i < Line.size(); i++)
			{
				generate_Hor_data[firstY][i] = Line[i];

			}
			if (initial_flag == 0)
			{
				if (outlier_flag.back() == 0)
				{
					initial_model = firstY;
					initial_flag = 1;
				}
			}
			//	toc = clock();
			//	printf("one Line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		}// vong lap while
		//toc = clock();
		//printf("All line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		toc = clock();
		printf("Finishe Intial fistLine Ver: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

		////5*****************************************Finished Initial first line*******************************************************************************************//
		//6*****************************************Draw Linear line model*******************************************************************************************//
				
		// tinh toan thoi gian cua Initial first line
			toc = clock();
		printf("Initial first line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		for (int firstY = initial_model; firstY < Hor_dataXY.size(); firstY++)// a= 105 = Hor_dataXY.size(), b=3138 =Ver_X[0].size()
		{
			//tic = clock();
			cout << "Interations:" << firstY << endl;
			indexX = 0;
			indexY = firstY;
			firstindex = Hor_dataXY[indexY][indexX];
			sliding_window = { (indexX + 1), indexY };
			vector<int> Line;
			Line.push_back(firstindex);  // Line ban dau =1, diem dau tien trong Hor_dataXY, sau do move sang diem moi

			i = 0;
			
			// check end line!
			while (1)
			{
				i = i + 1;
				if (i == Hor_dataXY[0].size())
				{
					break;
				}

				// Make a Line from 5 points by optimization projection errors
				if (!Line.empty())
				{
					index1 = Line.back();// -1so voi matlab, Vi line luu giu size cua Hor_XY, nhung ma chi so xuat phat tu 0 => -1
				}
				sliding_window = { indexX + i, indexY };// o day 6 0 thi matlab la 7 1
				index2 = Hor_dataXY[sliding_window[1]][sliding_window[0]];// nho hon matlab 1 so
				//cout<<"a"<<Hor_XY[index1][0];
				vector < double> v1 = Hor_XY[index1];// -0.157, 23.1686      23.1686  -0.093
				vector < double> v2 = Hor_XY[index2];// -0.093 23.1686

				if (have_predict == 1)
				{
					if (dist_P2P(predict_point, Hor_XY[index2]) < spacing_min)
					{
						Hor_XY[index2] = predict_point;
						if (Line.back() != index2)
						{
							Line.push_back(index2);
						}
						have_predict = 0;
						continue;

					}
					else
					{
						if (Hor_XY[index2][0] < predict_point.at(1))
						{
							Hor_XY[index2] = predict_point;
							if (Line.back() != index2)
							{
								Line.push_back(index2);
							}
							have_predict = 0;
							continue;
						}
						//clock_t tic, toc;
						//	tic = clock();
						Hor_XY.push_back(predict_point);// truong hop mat diem, thi se tao them 1 diem luu vao cuoi cua Hor_XY, va index cua no dc luu vao line, tu do de ve thi ve cac index chua trong line
						Line.push_back(Hor_XY[0].size());
						have_predict = 0;
						//			toc = clock();
						//		printf("Loading file aa: %f seconds", (double)(toc - tic) / CLOCKS_PER_SEC);
						continue;
					}

				}
				double d = dist_P2P(Hor_XY[index1], Hor_XY[index2]);// distance chinh xac
				//TESTCASE
				//cout << d; // 0.064
				if (d < spacing_min)
				{
					cout << "Warning : Constraints1(CT1): 0.04<distance" << endl;
					Hor_XY[index2] = Hor_XY[index1];
					if (((Hor_XY[index2][0]) + (Hor_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					continue;
				}
				else if (d > spacing_max) {
					cout << "Warning : Constraints1(CT1): distance<0.073" << endl;

					if (((Hor_XY[index2][0]) + (Hor_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					else
					{
						HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
						SetConsoleTextAttribute(h, FOREGROUND_RED | FOREGROUND_INTENSITY);
						cout << "Warning: outlier" << endl;
						//                     % generate the points and run the loop again
						//                  % can avoid two kinds of situation the outliers point
						SetConsoleTextAttribute(h, 15);
						pos_outlier = { i - 1, firstY };// phai tru 1 do xuat phat tu 0, ma i o day khoi tao la 1
						table_outlier.push_back(pos_outlier);
						/*clock_t tic, toc;*/
						//		tic = clock();
						Ref_model(predict_point, pos_outlier, Hor_dataXY, generate_Hor_data, Hor_XY, ConstraintsX);
						//predict_point = Ref_model(pos_outlier, Hor_dataXY, generate_Hor_data, Hor_XY, ConstraintsX);
						have_predict = 1;
						//toc = clock();
						//					printf("Ref model outside: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

						i = i - 1;
					}
					continue;

				}
				else if (Hor_XY[index2][0] < Hor_XY[index1][0])
				{
					cout << Hor_XY[index2][0];
					cout << Hor_XY[index1][0];
					cout << "Warning : Next point should in right side" << endl;
					continue;
				}

				if (Line.size() == 5)
				{
					vector<double>point_X(5);
					vector<double>point_Y(5);
					for (int i = 0; i < 5; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X[i] = Hor_XY[Line.at(i)][0];// test case: -0.157,-0.093,-0.0291,0.0349,0.0988
						point_Y[i] = Hor_XY[Line.at(i)][1];//test case: 5 so 23.1686

					}
					LL_model = Polyfit(point_X, point_Y, 1); // test case: 0 23.1686
				}
				else if (Line.size() > 5 && Line.size() <= 10)
				{
					if (dist_P2L(Hor_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						// Predict points//
						pos_outlier = { i - 1, firstY };
						table_outlier.push_back(pos_outlier);
						Ref_model(predict_point, pos_outlier, Hor_dataXY, generate_Hor_data, Hor_XY, ConstraintsX);
						have_predict = 1;
						i = i - 1;

						continue;
					};
					vector<double>point_X1;
					vector<double>point_Y1;  // moi vong lap deu ve lai

					for (int i = 0; i < Line.size(); i++) // Line tu 0 den 4, mat lab la 1 den 5
					{
						//Hor_XY[300];
						point_X1.push_back(Hor_XY[Line.at(i)][0]);
						point_Y1.push_back(Hor_XY[Line.at(i)][1]);
					}

					LL_model = Polyfit(point_X1, point_Y1, 1);

				}

				else if (Line.size() > 10)
				{
					if (dist_P2L(Hor_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						// Predict points//
						pos_outlier = { i - 1, firstY };
						table_outlier.push_back(pos_outlier);
						Ref_model(predict_point, pos_outlier, Hor_dataXY, generate_Hor_data, Hor_XY, ConstraintsX);
						have_predict = 1;
						i = i - 1;

						continue;
					}
					vector<double>point_X2(10);
					vector<double>point_Y2(10);
					for (int i = (Line.size() - 9); i < Line.size() + 1; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X2.at(i - (Line.size() - 9)) = Hor_XY[Line.at(i - 1)][0];
						point_Y2.at(i - (Line.size() - 9)) = Hor_XY[Line.at(i - 1)][1];
					}
					LL_model = Polyfit(point_X2, point_Y2, 1);
				}
				if (Line.back() != index2)
				{
					Line.push_back(index2);
				}

			}// vong lap while
			for (int i = 0; i < Line.size(); i++)
			{
				generate_Hor_data[firstY][i] = Line[i];

			}
			//toc = clock();
			//printf("one Line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		}// vong lap for
		/*toc = clock();
		printf("All line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);*/
		toc = clock();
		printf("Finishe Draw Linear line model: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		// *****************************************Finished Draw Linear line model*******************************************************************************************//
		// *****************************************Homography         *******************************************************************************************//
		
		vector<Point2d> src;
		vector<Point2d> HorXY_homo;
		for (int i = 0; i < Hor_XY.size(); i++)
		{
			//src.push_back({ Hor_XY[i].at(1), Hor_XY[i].at(0) });   // Chu y thu tu cua Hor_XY de test, dung khi la VerXy, can dao thu tu lai
			src.push_back({ Hor_XY[i].at(0), Hor_XY[i].at(1) }); // Horizoltal giu nguyen thu tu
		}

		//	perspectiveTransform(Hor_XY, HorXY_homo, H);
		perspectiveTransform(src, HorXY_homo, H);
		/*toc = clock();
		printf("Homography: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);*/
		toc = clock();
		printf("Finished Homography Hor: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		// *****************************************FinishedHomography *******************************************************************************************//

		//*********************************************Write output txt file*****************************************************************************//
		vector<Point2d> XYT_hom = HorXY_homo;
		Point2d hehe = HorXY_homo[20];
		//tic = clock();
		vector<vector<double> > generate_data(generate_Hor_data.size(), vector<double>(generate_Hor_data[0].size()));
		generate_data = generate_Hor_data;
		
		vector<vector<double>>bound_left(cross_point_LandR.size() / 2, vector<double>(cross_point_LandR[0].size()));
		for (int i = 0; i < cross_point_LandR.size() / 2; i++)
		{
			bound_left[i] = cross_point_LandR[i];//23
		}
		vector<vector<double>>bound_right((cross_point_LandR.size() - cross_point_LandR.size() / 2), vector<double>(cross_point_LandR[0].size()));
		for (int i = cross_point_LandR.size() / 2 ; i < cross_point_LandR.size(); i++)
		{
			bound_right[i - (cross_point_LandR.size() / 2 )] = cross_point_LandR[i];//177,24
		}
		ofstream writerXY("outputHor_XY.scr"); // mot vi du co mot so thu muc ko the tao file de ghi dc

		if (writerXY.is_open()) // buoc nay khi viet code rat quan trong, vi can phai kiem tra ghi dc hay o
		{
		/*	writerXY << "zoom" << endl;
			writerXY << "C" << endl;
			writerXY << "0,0" << endl;
			writerXY << "1" << endl;*/
			writerXY << "Pline" << endl;
			int temp = 0;
			for (int i = 0; i < generate_data.size(); i++)
			{

				for (int j = 0; j < generate_data[0].size(); j++)
				{
					if (generate_data[i][j] == 0&j>0)
					{
						break;
						
					}
					else if (src[generate_data[i][j]].x>bound_left[i][0] && src[generate_data[i][j]].x<bound_right[i][0])// toa do luu theo kieu [219  , 174] [219,175], luu diem theo y truoc x sau. Xet diem truoc homo ko phai diem homo
					{
						writerXY << XYT_hom[generate_data[i][j]].x << "," << XYT_hom[generate_data[i][j]].y << endl;
						
					}
					else
					{
						
						continue;// chay ra gia tri moi, neu dung break thi  no se thoat ra khoi vong for luon
					}

				}

				writerXY << "\nPline" << endl;
			}

			writerXY << endl;
			cout << "Write succesfully" << endl;
		}
		else
		{
			cout << "Error" << endl;
		}

		writerXY.close();
		toc = clock();
		printf("Finished Write output Hor trong loop: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		//*********************************************Finished Write output txt file*****************************************************************************//

	/*	toc = clock();
		printf("Writing file: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);*/
	}
	toc = clock();
	printf("Finished Hor ngoai loop: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	vector<vector<double>>temp(Ver_XY.size(), vector < double >(Ver_XY[0].size()));
	temp = Ver_XY;
	for (int i = 0; i < Ver_XY.size(); i++)
	{
	Ver_XY[i][0] = temp[i][1];
	Ver_XY[i][1] = temp[i][0];
	}

	if (switch_V_H == 0 || switch_V_H == 2)
	{

		/*vector<vector<double>>temp(Ver_XY.size(), vector < double >(Ver_XY[0].size()));
		temp = Ver_XY;
		for (int i = 0; i < Ver_XY.size(); i++)
		{
			Ver_XY[i][0] = temp[i][1];
			Ver_XY[i][1] = temp[i][0];
		}*/

		//5***********************************Process Ver select line ********************************************************************************//
		//6*****************************************Initial first line*******************************************************************************************//
		/*tic = clock();*/
		vector<vector<int>> table_outlier; // khoi tao = 1x1
		vector<int> outlier_flag;// = { 0 };
		int indexX;
		int indexY;
		int firstindex;
		int initial_model = 0;
		int initial_flag = 0;
		int firstY = -1;
		int have_predict = 0;
		int i;
		int index1;
		int index2;
		vector<double> predict_point(2);
		predict_point.reserve(2);
		vector<int>sliding_window(2);
		vector<int> pos_outlier(2);
		vector <double> LL_model = { 1, 2 };// linear line model tim bang polyfit sua sau
		vector<vector<double> > generate_Ver_data(Ver_dataXY.size(), vector<double>(Ver_dataXY[0].size()));
		while (1)
		{
			//tic = clock();
			if (initial_flag == 1)
			{
				firstY = firstY - 1;
			}
			else
			{
				firstY = firstY + 1;
			}
			if (firstY == -1)
			{
				break;
			}
			cout << "firstY =" << firstY << endl;
			indexX = 0;
			indexY = firstY;
			//outlier_flag[firstY]=0; // cai nay nen dung pushback
			outlier_flag.push_back(0);
			if (Ver_dataXY.size() == 0)
			{
				continue;
			}


			firstindex = Ver_dataXY[indexY][indexX];
			sliding_window = { (indexX + 1), indexY };
			vector<int> Line;
			Line.push_back(firstindex);  // Line ban dau =1, diem dau tien trong Ver_dataXY, sau do move sang diem moi

			i = 0;

			//	// check end line!
			while (1)
			{
				if (outlier_flag.back() == 1 && initial_flag == 0)
				{
					break;
				}
				i = i + 1;
				if (i == Ver_dataXY[0].size())
				{
					break;
				}

				//		// Make a Line from 5 points by optimization projection errors
				if (!Line.empty())
				{
					index1 = Line.back();// -1so voi matlab, Vi line luu giu size cua Ver_XY, nhung ma chi so xuat phat tu 0 => -1
				}
				sliding_window = { indexX + i, indexY };// o day 6 0 thi matlab la 7 1
				index2 = Ver_dataXY[sliding_window[1]][sliding_window[0]];// nho hon matlab 1 so

				if (have_predict == 1)
				{
					if (dist_P2P(predict_point, Ver_XY[index2]) < spacing_min)
					{
						Ver_XY[index2] = predict_point;
						if (Line.back() != index2)
						{
							Line.push_back(index2);
						}
						have_predict = 0;
						continue;

					}
					else
					{
						if (Ver_XY[index2][0] < predict_point.at(1))
						{
							Ver_XY[index2] = predict_point;
							if (Line.back() != index2)
							{
								Line.push_back(index2);
							}
							have_predict = 0;
							continue;
						}
						//clock_t tic, toc;
						//	tic = clock();
						Ver_XY.push_back(predict_point);// truong hop mat diem, thi se tao them 1 diem luu vao cuoi cua Ver_XY, va index cua no dc luu vao line, tu do de ve thi ve cac index chua trong line
						Line.push_back(Ver_XY[0].size());
						have_predict = 0;
						//			toc = clock();
						//		printf("Loading file aa: %f seconds", (double)(toc - tic) / CLOCKS_PER_SEC);
						continue;
					}

				}
				double d = dist_P2P(Ver_XY[index1], Ver_XY[index2]);// distance chinh xac
				//		//TESTCASE
				//		//cout << d; // 0.064
				if (d < spacing_min)
				{
					cout << "Warning : Constraints1(CT1): 0.04<distance" << endl;
					Ver_XY[index2] = Ver_XY[index1];
					if (((Ver_XY[index2][0]) + (Ver_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					continue;
				}
				else if (d > spacing_max) {
					cout << "Warning : Constraints1(CT1): distance<0.073" << endl;

					if (((Ver_XY[index2][0]) + (Ver_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					else
					{
						HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
						SetConsoleTextAttribute(h, FOREGROUND_RED | FOREGROUND_INTENSITY);
						cout << "Warning: outlier" << endl;
						//                     % generate the points and run the loop again
						//                  % can avoid two kinds of situation the outliers point
						SetConsoleTextAttribute(h, 15);
						//	outlier_flag(firstY) = 1;
						outlier_flag.back() = 1;
						//---------------------------------------------predict point----------------------------------------------------//
						if (initial_flag == 1)
						{


							pos_outlier = { i - 1, firstY };// phai tru 1 do xuat phat tu 0, ma i o day khoi tao la 1
							table_outlier.push_back(pos_outlier);

							InitialRef_model(predict_point, pos_outlier, Ver_dataXY, generate_Ver_data, Ver_XY, ConstraintsX);

							have_predict = 1;
							//toc = clock();
							//					printf("Ref model outside: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

							i = i - 1;
						}
						//---------------------------------------------predict point----------------------------------------------------//
					}
					continue;

				}
				else if (Ver_XY[index2][0] < Ver_XY[index1][0])
				{
					cout << Ver_XY[index2][0];
					cout << Ver_XY[index1][0];
					cout << "Warning : Next point should in right side" << endl;
					continue;
				}

				if (Line.size() == 5)
				{
					vector<double>point_X(5);
					vector<double>point_Y(5);
					for (int i = 0; i < 5; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X[i] = Ver_XY[Line.at(i)][0];// test case: -0.157,-0.093,-0.0291,0.0349,0.0988
						point_Y[i] = Ver_XY[Line.at(i)][1];//test case: 5 so 23.1686

					}
					LL_model = Polyfit(point_X, point_Y, 1); // test case: 0 23.1686
				}
				else if (Line.size() > 5 && Line.size() <= 10)
				{
					if (dist_P2L(Ver_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						//outlier_flag(firstY) = 1;
						outlier_flag.back() = 1;
						if (initial_flag == 1)
						{


							// Predict points//
							pos_outlier = { i - 1, firstY };
							table_outlier.push_back(pos_outlier);
							InitialRef_model(predict_point, pos_outlier, Ver_dataXY, generate_Ver_data, Ver_XY, ConstraintsX);
							have_predict = 1;
							i = i - 1;
						}
						continue;
					};
					vector<double>point_X1;
					vector<double>point_Y1;  // moi vong lap deu ve lai

					for (int i = 0; i < Line.size(); i++) // Line tu 0 den 4, mat lab la 1 den 5
					{
						//Ver_XY[300];
						point_X1.push_back(Ver_XY[Line.at(i)][0]);
						point_Y1.push_back(Ver_XY[Line.at(i)][1]);
					}

					LL_model = Polyfit(point_X1, point_Y1, 1);

				}

				else if (Line.size() > 10)
				{
					if (dist_P2L(Ver_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						//	outlier_flag(firstY) = 1;
						outlier_flag.back() = 1;
						if (initial_flag == 1)
						{


							// Predict points//
							pos_outlier = { i - 1, firstY };
							table_outlier.push_back(pos_outlier);
							InitialRef_model(predict_point, pos_outlier, Ver_dataXY, generate_Ver_data, Ver_XY, ConstraintsX);
							have_predict = 1;
							i = i - 1;
						}
						continue;
					}
					vector<double>point_X2(10);
					vector<double>point_Y2(10);
					for (int i = (Line.size() - 9); i < Line.size() + 1; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X2.at(i - (Line.size() - 9)) = Ver_XY[Line.at(i - 1)][0];
						point_Y2.at(i - (Line.size() - 9)) = Ver_XY[Line.at(i - 1)][1];
					}
					LL_model = Polyfit(point_X2, point_Y2, 1);
				}
				if (Line.back() != index2)
				{
					Line.push_back(index2);
				}

			}// vong lap while thu nhat chay tat ca cac diem tren 1 hang
			for (int i = 0; i < Line.size(); i++)
			{
				generate_Ver_data[firstY][i] = Line[i];

			}
			if (initial_flag == 0)
			{
				if (outlier_flag.back() == 0)// khong co outlier
				{
					initial_model = firstY;
					initial_flag = 1;
				}
			}
			//	toc = clock();
			//	printf("one Line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		}// vong lap while
		//toc = clock();
		//printf("All line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		toc = clock();
		printf("Finishe Intial fistLine Ver: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

		////5*****************************************Finished Initial first line*******************************************************************************************//
		//6*****************************************Draw Linear line model*******************************************************************************************//

		// tinh toan thoi gian cua Initial first line
		toc = clock();
		printf("Initial first line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		for (int firstY = initial_model; firstY < Ver_dataXY.size(); firstY++)// a= 105 = Ver_dataXY.size(), b=3138 =Ver_X[0].size()
		{
			//tic = clock();
			cout << "Interations:" << firstY << endl;
			indexX = 0;
			indexY = firstY;
			firstindex = Ver_dataXY[indexY][indexX];
			sliding_window = { (indexX + 1), indexY };
			vector<int> Line;
			Line.push_back(firstindex);  // Line ban dau =1, diem dau tien trong Ver_dataXY, sau do move sang diem moi

			i = 0;

			// check end line!
			while (1)
			{
				i = i + 1;
				if (i == Ver_dataXY[0].size())
				{
					break;
				}

				// Make a Line from 5 points by optimization projection errors
				if (!Line.empty())
				{
					index1 = Line.back();// -1so voi matlab, Vi line luu giu size cua Ver_XY, nhung ma chi so xuat phat tu 0 => -1
				}
				sliding_window = { indexX + i, indexY };// o day 6 0 thi matlab la 7 1
				index2 = Ver_dataXY[sliding_window[1]][sliding_window[0]];// nho hon matlab 1 so
				//cout<<"a"<<Ver_XY[index1][0];
				vector < double> v1 = Ver_XY[index1];// -0.157, 23.1686      23.1686  -0.093
				vector < double> v2 = Ver_XY[index2];// -0.093 23.1686

				if (have_predict == 1)
				{
					if (dist_P2P(predict_point, Ver_XY[index2]) < spacing_min)
					{
						Ver_XY[index2] = predict_point;
						if (Line.back() != index2)
						{
							Line.push_back(index2);
						}
						have_predict = 0;
						continue;

					}
					else
					{
						if (Ver_XY[index2][0] < predict_point.at(1))
						{
							Ver_XY[index2] = predict_point;
							if (Line.back() != index2)
							{
								Line.push_back(index2);
							}
							have_predict = 0;
							continue;
						}
						//clock_t tic, toc;
						//	tic = clock();
						Ver_XY.push_back(predict_point);// truong hop mat diem, thi se tao them 1 diem luu vao cuoi cua Ver_XY, va index cua no dc luu vao line, tu do de ve thi ve cac index chua trong line
						Line.push_back(Ver_XY[0].size());
						have_predict = 0;
						//			toc = clock();
						//		printf("Loading file aa: %f seconds", (double)(toc - tic) / CLOCKS_PER_SEC);
						continue;
					}

				}
				double d = dist_P2P(Ver_XY[index1], Ver_XY[index2]);// distance chinh xac
				//TESTCASE
				//cout << d; // 0.064
				if (d < spacing_min)
				{
					cout << "Warning : Constraints1(CT1): 0.04<distance" << endl;
					Ver_XY[index2] = Ver_XY[index1];
					if (((Ver_XY[index2][0]) + (Ver_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					continue;
				}
				else if (d > spacing_max) {
					cout << "Warning : Constraints1(CT1): distance<0.073" << endl;

					if (((Ver_XY[index2][0]) + (Ver_XY[index2][1])) == 0)
					{
						cout << "Warning : XY(index2,:)==0" << endl;
						break;
					}
					else
					{
						HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
						SetConsoleTextAttribute(h, FOREGROUND_RED | FOREGROUND_INTENSITY);
						cout << "Warning: outlier" << endl;
						//                     % generate the points and run the loop again
						//                  % can avoid two kinds of situation the outliers point
						SetConsoleTextAttribute(h, 15);
						pos_outlier = { i - 1, firstY };// phai tru 1 do xuat phat tu 0, ma i o day khoi tao la 1
						table_outlier.push_back(pos_outlier);
						/*clock_t tic, toc;*/
						//		tic = clock();
						Ref_model(predict_point, pos_outlier, Ver_dataXY, generate_Ver_data, Ver_XY, ConstraintsX);
						//predict_point = Ref_model(pos_outlier, Ver_dataXY, generate_Ver_data, Ver_XY, ConstraintsX);
						have_predict = 1;
						//toc = clock();
						//					printf("Ref model outside: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

						i = i - 1;
					}
					continue;

				}
				else if (Ver_XY[index2][0] < Ver_XY[index1][0])
				{
					cout << Ver_XY[index2][0];
					cout << Ver_XY[index1][0];
					cout << "Warning : Next point should in right side" << endl;
					continue;
				}

				if (Line.size() == 5)
				{
					vector<double>point_X(5);
					vector<double>point_Y(5);
					for (int i = 0; i < 5; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X[i] = Ver_XY[Line.at(i)][0];// test case: -0.157,-0.093,-0.0291,0.0349,0.0988
						point_Y[i] = Ver_XY[Line.at(i)][1];//test case: 5 so 23.1686

					}
					LL_model = Polyfit(point_X, point_Y, 1); // test case: 0 23.1686
				}
				else if (Line.size() > 5 && Line.size() <= 10)
				{
					if (dist_P2L(Ver_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						// Predict points//
						pos_outlier = { i - 1, firstY };
						table_outlier.push_back(pos_outlier);
						Ref_model(predict_point, pos_outlier, Ver_dataXY, generate_Ver_data, Ver_XY, ConstraintsX);
						have_predict = 1;
						i = i - 1;

						continue;
					};
					vector<double>point_X1;
					vector<double>point_Y1;  // moi vong lap deu ve lai

					for (int i = 0; i < Line.size(); i++) // Line tu 0 den 4, mat lab la 1 den 5
					{
						//Ver_XY[300];
						point_X1.push_back(Ver_XY[Line.at(i)][0]);
						point_Y1.push_back(Ver_XY[Line.at(i)][1]);
					}

					LL_model = Polyfit(point_X1, point_Y1, 1);

				}

				else if (Line.size() > 10)
				{
					if (dist_P2L(Ver_XY[index2], LL_model) > dist_LL_Model)
					{
						cout << "Warning : Base on current linear line model, projection distance < Threshold" << endl;
						// Predict points//
						pos_outlier = { i - 1, firstY };
						table_outlier.push_back(pos_outlier);
						Ref_model(predict_point, pos_outlier, Ver_dataXY, generate_Ver_data, Ver_XY, ConstraintsX);
						have_predict = 1;
						i = i - 1;

						continue;
					}
					vector<double>point_X2(10);
					vector<double>point_Y2(10);
					for (int i = (Line.size() - 9); i < Line.size() + 1; i++) // Line tu 0 den 4, mat lab la 1 den 5
					{

						point_X2.at(i - (Line.size() - 9)) = Ver_XY[Line.at(i - 1)][0];
						point_Y2.at(i - (Line.size() - 9)) = Ver_XY[Line.at(i - 1)][1];
					}
					LL_model = Polyfit(point_X2, point_Y2, 1);
				}
				if (Line.back() != index2)
				{
					Line.push_back(index2);
				}

			}// vong lap while
			for (int i = 0; i < Line.size(); i++)
			{
				generate_Ver_data[firstY][i] = Line[i];

			}
			//toc = clock();
			//printf("one Line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		}// vong lap for
		/*toc = clock();
		printf("All line: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);*/
		toc = clock();
		printf("Finishe Draw Linear line model: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		// *****************************************Finished Draw Linear line model*******************************************************************************************//
		// *****************************************Homography         *******************************************************************************************//
		//	vector<vector<double> > VerXY_homo(Ver_XY.size(), vector<double>(Ver_XY[0].size()));
		vector<Point2d> src;
		vector<Point2d> VerXY_homo;
		for (int i = 0; i < Ver_XY.size(); i++)
		{
			src.push_back({ Ver_XY[i].at(1), Ver_XY[i].at(0) });   // Chu y thu tu cua Ver_XY de test, dung khi la VerXy, can dao thu tu lai
			//	src.push_back({ Ver_XY[i].at(0), Ver_XY[i].at(1) }); // Horizoltal giu nguyen thu tu
		}
		//vector<Point2d> src;
		//vector<Point2d> HorXY_homo;
		//for (int i = 0; i < Hor_XY.size(); i++)
		//{
		//	//src.push_back({ Hor_XY[i].at(1), Hor_XY[i].at(0) });   // Chu y thu tu cua Hor_XY de test, dung khi la VerXy, can dao thu tu lai
		//	src.push_back({ Hor_XY[i].at(0), Hor_XY[i].at(1) }); // Horizoltal giu nguyen thu tu
		//}
		//	perspectiveTransform(Ver_XY, VerXY_homo, H);
		perspectiveTransform(src, VerXY_homo, H);
		/*toc = clock();
		printf("Homography: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);*/
		toc = clock();
		printf("Finish Homorgraph Ver: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		// *****************************************FinishedHomography *******************************************************************************************//

		//*********************************************Write output txt file*****************************************************************************//
		vector<Point2d> XYT_hom = VerXY_homo;
		Point2d hehe = VerXY_homo[20];
		//tic = clock();
		vector<vector<double> > generate_data(generate_Ver_data.size(), vector<double>(generate_Ver_data[0].size()));
		generate_data = generate_Ver_data;


		vector<vector<double>>bound_down(cross_point_UandD.size() / 2, vector<double>(cross_point_UandD[0].size()));
		for (int i = 0; i < cross_point_UandD.size() / 2; i++)
		{
			bound_down[i] = cross_point_UandD[i];
		}
		vector<vector<double>>bound_up((cross_point_UandD.size() - cross_point_UandD.size() / 2), vector<double>(cross_point_UandD[0].size()));
		for (int i = cross_point_UandD.size() / 2 ; i < cross_point_UandD.size(); i++)
		{
			bound_up[i - (cross_point_UandD.size() / 2 )] = cross_point_UandD[i];
		}

		ofstream writerXY("outputVER_XY.scr"); // mot vi du co mot so thu muc ko the tao file de ghi dc
		int ee = 0;
		if (writerXY.is_open()) // buoc nay khi viet code rat quan trong, vi can phai kiem tra ghi dc hay o
		{
			/*writerXY << "zoom" << endl;
			writerXY << "C" << endl;
			writerXY << "0,0" << endl;
			writerXY << "1" << endl;*/
			writerXY << "Pline" << endl;
			int temp = 0;
			for (int i = 0; i < generate_data.size(); i++)
			{
				//	writerXY << "\nPline" << endl;
				for (int j = 0; j < generate_data[0].size(); j++)
				{
					if (generate_data[i][j] == 0 & j>0)// do co truong hop dau tien j=0
					{
						break;
						//	writerXY << XYT_hom[generate_data[i][j]].x << "," << XYT_hom[generate_data[i][j]].y << endl;
					}
					else if (src[generate_data[i][j]].y>bound_down[i][1] && src[generate_data[i][j]].y<bound_up[i][1])
					{
						writerXY << XYT_hom[generate_data[i][j]].x << "," << XYT_hom[generate_data[i][j]].y << endl;
					}
					else
					{
						//	writerXY << XYT_hom[generate_data[i][j]].x << "," << XYT_hom[generate_data[i][j]].y << endl;
						//cout <<"\na:"<< src[generate_data[i][j]].y;
						continue;
					}

				}

				writerXY << "\nPline" << endl;
			}

			writerXY << endl;
			cout << "Write succesfully" << endl;
		}
		else
		{
			cout << "Error" << endl;
		}

		writerXY.close();
		toc = clock();
		printf("Finish write output Ver in boundary infile trong loop: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		//*********************************************Finished Write output txt file*****************************************************************************//

		/*	toc = clock();
		printf("Writing file: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);*/
	}
	toc = clock();
	printf("Finished Ver ngoai loop: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

	
	
	StartProgram1("C:\\Program Files\\Autodesk\\AutoCAD 2012 - English\\acad.exe\ ", "/b  D:\\Project\\Opencv20_5_Full\\HelloWorld\\outputHor_XY.scr ");// location of the file
	StartProgram1("C:\\Program Files\\Autodesk\\AutoCAD 2012 - English\\acad.exe\ ", "/b  D:\\Project\\Opencv20_5_Full\\HelloWorld\\outputVER_XY.scr ");//location of the file
	
		toc = clock();
		cout << "\n tic:"<<tic;
		printf("All: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		system("pause");
		
		return 0;
	}


	