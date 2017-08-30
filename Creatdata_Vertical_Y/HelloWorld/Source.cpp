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


using namespace std;
//using std::vector;
const double spacing_min = 0.04;
const double spacing_max = 0.073;
const double ConstraintsX = 0.065;
const double dist_LL_Model = 0.03;
int switch_V_H = 2;// 0 : vertical, 1 : horizol, 2: both

void StartProgram1(string filename, string filename1)
{
	ShellExecuteA(NULL, "open", filename.c_str(), filename1.c_str(), NULL, SW_SHOW);

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
				printf("hang:  %d", countLines);
				printf("so cot:  %d\n", columns);
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

	int expand = 150;
	clock_t tic, toc;
	tic = clock();
	cout << "tic:" << tic << endl;

	int a, b;
	int a_H, b_H;
	int *p;
	//**************************Read data******************************************************************************************************//
	toc = clock();
	printf("time to homography: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	//// Vertical data
	p = readFile("105L.V.G2X.TXT", &a, &b);//105L.V.G1X.TXT
	vector<vector<double> > Ver_X(a, vector<double>(b));
	toc = clock();
	printf("read number of rows,columns of Vertical: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	
	// Luu gia tri Y*****************************************************************************//
	vector<vector<double> > Ver_Y(a, vector<double>(b));
	// Đọc các giá trị của ma trận từ file
	ifstream V_L("105L.V.G2Y.TXT");//105L.V.G1Y.TXT

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
	// ********** write newfile Vertical

	ofstream writerXYZ("testdata105L.V.G2Y_150.txt"); // mot vi du co mot so thu muc ko the tao file de ghi dc
	//double temp = 0.064 * 23;
	//double temp = Ver_Y[104][0] - Ver_Y[0][0] + Ver_Y[1][0] - Ver_Y[0][0];
	if (writerXYZ.is_open()) // buoc nay khi viet code rat quan trong, vi can phai kiem tra ghi dc hay o
	{
		for (int i = 0; i< (Ver_Y.size()+expand); i++)//ee so muon mo rong theo phuong ngang
		{
			for (int j = 0; j < Ver_Y[0].size() ; j++)
			{

				if (i<=104)// giu nguyen cac gia tri sau 2750
				{
					if (j == Ver_Y[0].size()-1)
					{
						writerXYZ << Ver_Y[i][j];
						continue;
					}
					else
					{
						writerXYZ << Ver_Y[i][j] << " ";
						continue;
					}
					
				}
				else// m>104 thi bang 104
				{
					if (j == Ver_Y[0].size() - 1)
					{
						writerXYZ << Ver_Y[104][j];
						continue;
					}
					else
					{
						writerXYZ << Ver_Y[104][j] << " ";
						continue;
					}
					
				}
				
			}
			writerXYZ << endl;
		}		

	}
	else
	{
		cout << "Error" << endl;
	}


	writerXYZ.close();
	// finished write new file

	// finished write new file
	// read new data
	int testa, testb;

	p = readFile("testdata105L.V.G2Y_150.txt", &testa, &testb);//105L.H.G1X.TXT
	vector<vector<double> > testVer_Y(testa, vector<double>(testb));
	// Đọc các giá trị của ma trận từ file
	ifstream testH_V_L("testdata105L.V.G2Y_150.TXT");//105L.H.G1Y.TXT

	if (testH_V_L.is_open())
	{

		for (int i = 0; i < testa; i++)// a=105 b= 3138
		{
			for (int j = 0; j < testb; j++)
			{
				// xoa nhung gia tri lon hon so column cua dong do
				if (j> *(p + i))
				{

					testVer_Y[i][j] = 0;
					continue;
				}
				testH_V_L >> testVer_Y[i][j];
			}


		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}
	testH_V_L.close();
	toc = clock();
	printf("Store Horizoltal_Y data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	

	//1************************** finished Read data******************************************************************************************************//


	StartProgram1("C:\\Program Files\\Autodesk\\AutoCAD 2012 - English\\acad.exe\ ", "/b  D:\\Project\\Opencv20_5_Full\\HelloWorld\\outputHor_XY.scr ");
	StartProgram1("C:\\Program Files\\Autodesk\\AutoCAD 2012 - English\\acad.exe\ ", "/b  D:\\Project\\Opencv20_5_Full\\HelloWorld\\outputVER_XY.scr ");

	toc = clock();
	cout << "\n tic:" << tic;
	printf("All: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	system("pause");

	return 0;
}


