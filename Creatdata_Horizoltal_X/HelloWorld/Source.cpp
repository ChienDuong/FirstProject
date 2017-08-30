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


	clock_t tic, toc;
	tic = clock();
	cout << "tic:" << tic << endl;

	int a, b;
	int a_H, b_H;
	int *p;
	//**************************Read data******************************************************************************************************//
	toc = clock();
	int expand = 150;// he so mo rong
	printf("time to homography: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	// Horizoltal data
	p = readFile("105L.H.G2X.TXT", &a_H, &b_H);//105L.H.G1X.TXT
	vector<vector<double> > Hor_X(a_H, vector<double>(b_H));
	toc = clock();
	printf("read number of rows,columns: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	// Đọc các giá trị của ma trận từ file
	ifstream H_H_L("105L.H.G2X.TXT");//105L.H.G1X.TXT

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

	cout << "Read succesfully" << endl;

	H_H_L.close();
	toc = clock();
	printf("Store Horizoltal_X data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

	// write newfile

	ofstream writerXY("testdata105L.H.G2X_150.txt"); // mot vi du co mot so thu muc ko the tao file de ghi dc
	double temp = 0.064;
	if (writerXY.is_open()) // buoc nay khi viet code rat quan trong, vi can phai kiem tra ghi dc hay o
	{

		for (int m = 0; m < Hor_X.size(); m++)//ee so muon mo rong theo phuong ngang
		{
			for (int j = 0; j < Hor_X[0].size() + 23*expand; j++)
			{

				if (j == Hor_X[0].size() + 23*expand - 1)
				{
					if (Hor_X[m][2750 + (j - 2750 - 23*expand)] == 0)
					{
						writerXY << 0 ;
						break;
					}
					else
					{
						writerXY << Hor_X[m][2750 + (j - 2750 - 23*expand)] + temp * 23*expand;
						break;
					}
					
					
				}
				if (j <= 2750)
				{
					writerXY << Hor_X[m][j] << " ";
					continue;
				}
				if (j > 2750 && j <= 2750 + 23*expand)
				{
					writerXY << Hor_X[m][2750] + temp*(j - 2750) << " ";
					continue;
				}
				if (j > (2750 + 23*expand))
				{
					if (Hor_X[m][2750 + (j - 2750 - 23*expand)]==0)
					{
						writerXY << 0 << " ";
						continue;
					}
					else
					{
						writerXY << Hor_X[m][2750 + (j - 2750 - 23*expand)] + temp * 23*expand << " ";
						continue;
					}
					
				}
				else
				{
					cout << "Error!";
					getchar();
				}

				//writerXY << Hor_X[i - Hor_X.size()*k][j - *(p + (i - Hor_X.size()*k))*m] + temp*k << " ";// k la so lap lai hang, m la so lap lai cot


			}
			writerXY << endl;
		}

	}
	else
	{
		cout << "Error" << endl;
	}


	writerXY.close();
	// finished write new file



	// read new data
	int testa_H, testb_H;

	p = readFile("testdata105L.H.G2X_150.txt", &testa_H, &testb_H);//105L.H.G1X.TXT
	vector<vector<double> > testHor_X(testa_H, vector<double>(testb_H));
	// Đọc các giá trị của ma trận từ file
	ifstream testH_H_L("testdata105L.H.G2X_150.TXT");//105L.H.G1Y.TXT

	if (testH_H_L.is_open())
	{

		for (int i = 0; i < testa_H; i++)// a=105 b= 3138
		{
			for (int j = 0; j < testb_H; j++)
			{
				// xoa nhung gia tri lon hon so column cua dong do
				if (j> *(p + i))
				{

					testHor_X[i][j] = 0;
					continue;
				}
				testH_H_L >> testHor_X[i][j];
			}


		}
	}
	else
	{
		cout << "Error! Cannot open file!" << endl;
	}
	testH_H_L.close();
	toc = clock();
	printf("Store Horizoltal_Y data in: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

	
	////finished readnew data
	
	//1************************** finished Read data******************************************************************************************************//


	StartProgram1("C:\\Program Files\\Autodesk\\AutoCAD 2012 - English\\acad.exe\ ", "/b  D:\\Project\\Opencv20_5_Full\\HelloWorld\\outputHor_XY.scr ");
	StartProgram1("C:\\Program Files\\Autodesk\\AutoCAD 2012 - English\\acad.exe\ ", "/b  D:\\Project\\Opencv20_5_Full\\HelloWorld\\outputVER_XY.scr ");

	toc = clock();
	cout << "\n tic:" << tic;
	printf("All: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	system("pause");

	return 0;
}


