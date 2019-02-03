// Lab1Slau.cpp: ���������� ����� ����� ��� ����������� ����������.
//

#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <malloc.h>
#include <cmath>
using namespace std;
typedef double mytype;

mytype** readfromfile(string filename, int &n)
{
	char buff[500];													
	mytype element;
	int count = 0;
	ifstream fin(filename); // ������� ���� ��� ������		

	while (!fin.eof())//������������ ����� �����
	{
		fin.getline(buff, 500);
		count++;
	}
	//�������� ������ ��� ������
	int i;
	mytype **extMatrix;
	extMatrix = new mytype*[count];
	for (i = 0; i<count; i++)
		//a[i], a[i] �������� � ��������� ���� double
		extMatrix[i] = new mytype[count + 1];

	fin.close(); // ��������� ����
	fin.open(filename); // ������� ���� ��� ������
	int elcnt = 0;
	while (!fin.eof())
	{
		fin >> element; // ������� ��������� �������
		div_t ij = div(elcnt, count + 1);
		extMatrix[ij.quot][ij.rem] = element;
		elcnt++;
	}
	n = count;
	fin.close(); // ��������� ����
	return extMatrix;
}

void matrix_destroyer(mytype** ary, int n)
{
	for (int i = 0; i < n; i++) {
		delete[] ary[i];
	}
	delete[] ary;
}

mytype** readfromscreen(int &k)
{
	int i, j, raz;
	cout << "������� ����������� �������: ";
	cin >> raz;
	k = raz;
	mytype **extMatrix;
	extMatrix = new mytype*[raz];

	for (i = 0; i<raz; i++)
		extMatrix[i] = new mytype[raz + 1];

	for (i = 0; i < raz; i++)
	{
		for (j = 0; j < raz+1; j++)
		{
 			cin >> extMatrix[i][j];
		}
	}
	return extMatrix;
}

mytype** readfromscreensquare(int &k)
{
	int i, j, raz;
	cout << "������� ����������� ������� A: ";
	cin >> raz;
	k = raz;
	mytype **extMatrix;
	extMatrix = new mytype*[raz];

	for (i = 0; i<raz; i++)
		//a[i], a[i] �������� � ��������� ���� double
		extMatrix[i] = new mytype[raz];

	for (i = 0; i < raz; i++)
	{
		for (j = 0; j < raz; j++)
		{
			cin >> extMatrix[i][j];
		}
	}
	return extMatrix;
}

void printMatrix(mytype** extMatrix, int k, int m, bool extended)
{
	for ( int i = 0; i < k; i++)
	{
		for (int j = 0; j < m; j++)
		{
			char simb;
			if ((j == k - 1) && extended) { simb = '='; }
			else { simb = ' '; }
			cout << extMatrix[i][j]<< simb;
	
		}
		cout << endl;
	}

}

mytype** getMatrixMinor(mytype **extMatrix, int i, int j, int n)
{
	mytype **minor;
	minor = new mytype*[n - 1];

	for (int ii = 0; ii<n - 1; ii++)
		//a[i], a[i] �������� � ��������� ���� double
		minor[ii] = new mytype[n - 1];

	int ki, kj, di, dj;
	di = 0;
	for (ki = 0; ki< n - 1; ki++)
	{ // �������� ������� ������
		if (ki == i) di = 1;
		dj = 0;

		for (kj = 0; kj< n - 1; kj++)
		{ // �������� ������� �������
			if (kj == j) dj = 1;
			minor[ki][kj] = extMatrix[ki + di][kj + dj];
		}
	}
	return minor;
}

mytype determinant(mytype **extMatrix, int m)
{
	int i, j, k, n;
	mytype det = 0;
	mytype **minor;
	minor = new mytype*[m-1];

	for (i = 0; i<m; i++)
		minor[i] = new mytype[m-1];
	j = 0;
	k = 1; //(-1) � ������� i
	n = m - 1;

	if (m<1) cout << "������������ ��������� ����������!";

	if (m == 1)
	{
		det = extMatrix[0][0];
		return det;
	}

	if (m>=2)
	{
		for (i = 0; i<m; i++) {
			minor = getMatrixMinor(extMatrix, i, 0, m);
			//cout << extMatrix[i][j] << endl;
			//PrintMatr(p, n);
			det = det + k * extMatrix[i][0] * determinant(minor, n);
			k = -k;
		}
	}
	matrix_destroyer(minor, n);//������� ������ ���������� ��� �����
	return det;
}

mytype** multiplyMatrix(mytype** Matrix1, mytype** Matrix2, int n)
{
	mytype **result;
	result = new mytype*[n];
	mytype s = 0;
	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		for (int l = 0; l < n; l++)
		{
			s = 0;

			for (int j = 0; j < n; j++)
			{
				s += Matrix1[i][j] * Matrix2[j][l];
			}
			result[i][l] = s;
		}
	}
	return result;
}

mytype* multiplyMatrixVector(mytype** Matrix, mytype* Vector, int n)
{
	mytype *result;
	result = new mytype[n];
	mytype s = 0;
	
	for (int i = 0; i < n; i++)
	{
		s = 0;

		for (int j = 0; j < n; j++)
		{
			s += Matrix[i][j] * Vector[j];
		}
		result[i] = s;
		
	}
	return result;
}

mytype** transposeMatrix(mytype** Matrix, int n)
{
	mytype **result;
	result = new mytype*[n];
	mytype boof;
	for (int ii = 0; ii<n; ii++)
	result[ii] = new mytype[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			result[i][j] = Matrix[j][i];
			result[j][i] = Matrix[i][j];
		}
	}

	return result;
}

void unitMatrix(mytype** Matrix, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j) { Matrix[i][j] = 1; }
			else  Matrix[i][j] = 0 ;
		}
	}
}

mytype* gaussmethod(mytype** extMatrix, int n)
{
	mytype *solution;
	solution = new mytype[n];
	mytype maxvalue = 0;
	int imax;

	for (int cnt = 0; cnt < n; cnt++)
	{
		solution[cnt] = 0;
	}
	/*double det = determinant(extMatrix, n);
	if (abs(det) < 1e-30) 
	{
		cout << "������������ ����� 0. �� ���������� ������������� �������." << endl;
		solution = nullptr;
	}
	else
	{*/

	for (int i = 0; i < n - 1; i++)//���� �� �������, ������� ���������� �� �����������
	{
		//����� ���� �������� �� i-�� �������
		maxvalue = 0;
		for (int il = i; il < n; il++)
		{
			if (maxvalue < abs(extMatrix[il][i]))
			{
				maxvalue = abs(extMatrix[il][i]);
				imax = il;
			}
		}

		if (maxvalue < 1e-10)
		{
			cout << "�� ���������� ������������� �������." << endl;
			return nullptr;
		}

		if (imax != i)
		{
			mytype* buf = extMatrix[imax];
			extMatrix[imax] = extMatrix[i];
			extMatrix[i] = buf;
		}

		//extMatrix[i][n] = extMatrix[i][n] / extMatrix[i][i];
		mytype aii = extMatrix[i][i];

		if (abs(aii) < 1e-10)
		{
			cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
			return nullptr;
		}
	
			for (int j = i; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
			{
				extMatrix[i][j] = extMatrix[i][j] / aii;
			}

			for (int ii = i + 1; ii < n; ii++)//��������� �� ���������� ����� i-�� ������
			{
				mytype a_ii_i = extMatrix[ii][i];
				for (int jj = i; jj <= n; jj++)
				{
					extMatrix[ii][jj] -= a_ii_i * extMatrix[i][jj];
				}
			}
	}
			//��������� ������ ������
	mytype	 aii = extMatrix[n - 1][n - 1];
			if (abs(aii) < 1e-10)
			{
				cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
				return nullptr;
			}
			for (int j = n - 1; j <= n; j++)//���� �� ��������� ��������, ������� ���������� �� �����������  �� i+1???
			{
				extMatrix[n - 1][j] = extMatrix[n - 1][j] / aii;
			}
			//printMatrix(extMatrix, n, n + 1, true);
			//�������� ���

			mytype sum = 0;
			for (int i = n - 1; i >= 0; i--)
			{
				sum = 0;
				for (int j = i + 1; j < n; j++) //��������� ��� ����� ������� ����������  ���������� �� ������������ ������� ������
				{
					sum += solution[j] * extMatrix[i][j];
				}
				solution[i] = extMatrix[i][n] - sum;//�������� �� ������ ����� 
			}
	
	//printMatrix(extMatrix, n);//������ ������������������� (��� ��������)
	return solution;
}
 
mytype* qrmethod(mytype** extMatrix, int n, bool flag)
{
	mytype c, s;
	mytype *solution;
	solution = new mytype[n];
	mytype **t;
	t = new mytype*[n];
	mytype **tij;
	tij = new mytype*[n];
	mytype **r;
	r = new mytype*[n];
	mytype *b;
	b = new mytype[n];
	mytype **eiimax;
	eiimax = new mytype*[n];
	mytype **multiplyqr;
	multiplyqr = new mytype*[n];

	mytype **oldmemory;

	for (int ii = 0; ii < n; ii++)
	{
		t[ii] = new mytype[n];
		tij[ii] = new mytype[n];
		r[ii] = new mytype[n];
		eiimax[ii] = new mytype[n];
		multiplyqr[ii] = new mytype[n];
	}

	for (int i = 0; i < n; i++)//����������� ��������� ���� �� ������� � � ������ b
	{
		b[i] = extMatrix[i][n];
	}

	for (int i = 0; i < n; i++)//��������� �������� �� ������� ������� � ������� r
	{
		for (int j = 0; j < n; j++)
		{		
			r[i][j] = extMatrix[i][j];
		}
	}

	unitMatrix(t, n);

	mytype maxvalue = 0;
	int imax;

	for (int i = 0; i < n - 1; i++)
	{
	
		maxvalue = 0;
		for (int il = i; il < n; il++)
		{
			if (maxvalue < abs(extMatrix[il][i]))
			{
				maxvalue = abs(extMatrix[il][i]);
				imax = il;
			}
		}

		if (maxvalue < 1e-10)
		{
			cout << "�� ���������� ������������� �������." << endl;
			return nullptr;
		}

		if (imax != i)
		{//������������ �����, � i-�� ������ ���������� imax-������ � ��������
			mytype* buf = r[imax];
			r[imax] = r[i];
			r[i] = buf;
        //��������� ������� t �� ��������������� ������� eiimax

			unitMatrix(eiimax, n);//�������������� ������� t ��� ��������� ��� ������������ ������������
		
			eiimax[i][i] = 0;//��������� ��������������� �������
			eiimax[imax][imax] = 0;
			eiimax[i][imax] = 1;
			eiimax[imax][i] = 1;

			oldmemory = t;
			t = multiplyMatrix(eiimax, t, n);
			matrix_destroyer(oldmemory, n);
		}
		
		
		for (int j = i + 1; j < n; j++)
		{
			c = r[i][i] / (sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]));
			s = r[j][i] / (sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]));
		
			unitMatrix(tij, n);//tij ���������������� ��������� ��������

			tij[i][i] = c;
			tij[i][j] = s;
			tij[j][i] = -s;
			tij[j][j] = c;
			/*cout << "������� t" <<i<<j<< endl;
			printMatrix(tij, n, n, false);*/
			oldmemory = t;
			t = multiplyMatrix(tij, t, n);//��������� ������� t �� ������ ����� ��������
			matrix_destroyer(oldmemory, n);
			oldmemory = r;
			r = multiplyMatrix(tij, r, n);//��������� ������� r �� ������ ����� ��������
			matrix_destroyer(oldmemory, n);
			/*cout << "������� r" << i << j << endl;
			printMatrix(r, n, n, false);*/			
		}
	}

	//��������� ������ ������ �� ��� ������� ������������
	mytype	 aii = r[n - 1][n - 1];
	if (abs(aii) < 1e-10)
	{
		cout << "�� ���������� ������������� �������. ��������� ������ ������������������� ������� - �������" << endl;
		return nullptr;
	}

	mytype** q = transposeMatrix(t, n);
	if (flag)
	{
		cout << endl;
		cout << "������� Q:" << endl;
		cout << endl;
		printMatrix(q, n, n, false);
		cout << endl;
		cout << "������� R:" << endl;
		cout << endl;
		printMatrix(r, n, n, false);
		cout << endl;
		cout << "������������ ������ Q � R:" << endl;
		cout << endl;
		multiplyqr = multiplyMatrix(q, r, n);
		printMatrix(multiplyqr, n, n, false);
	}

	b = multiplyMatrixVector(t, b, n);//����������� ������ ��������� ������
	for (int cnt = 0; cnt < n; cnt++)
	{
		solution[cnt] = 0;
	}

	mytype sum = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		sum = 0;
		for (int j = i + 1; j < n; j++) //��������� ��� ����� ������� ����������  ���������� �� ������������ ������� ������
		{
			sum += solution[j] * r[i][j];
		}
		solution[i] = (b[i] - sum)/ r[i][i];//�������� �� ������ ����� 
	}

	return solution;
}

mytype discrepancy(mytype** extMatrix, mytype* solution, int n)
{
	mytype *b1;
	b1 = new mytype[n];
	mytype result = 0;

	b1 = multiplyMatrixVector(extMatrix, solution, n);
	
	for (int i = 0; i < n-1; i++)
	{
		result += pow((b1[i]-extMatrix[i][n]), 2);
	}
	result = sqrt(result);
	return result;
}

mytype** reverseMatrix(mytype** Matrix, int n)
{
	mytype **extMatrix;
	extMatrix = new mytype*[n];
	mytype **result;
	result = new mytype*[n];
	mytype *ordinary;
	ordinary = new mytype[n];

	for (int ii = 0; ii<n; ii++)
		result[ii] = new mytype[n];
	for (int ii = 0; ii<n; ii++)
	    extMatrix[ii] = new mytype[n+1];

		for (int j = 0; j < n; j++)//��������� ������ ������� extMatrix ������
		{
			extMatrix[j][n] = 0;
		}

	for (int i = 0; i < n; i++)
	{
		for (int ii = 0; ii < n; ii++)//��������� �������� �� ������� � ����������� �������
		{
			for (int jj = 0; jj < n; jj++)
			{
				extMatrix[ii][jj] = Matrix[ii][jj];
			}
		}
		for (int jj = 0; jj < n; jj++)//��������� ������ ������� extMatrix ������
		{
			extMatrix[jj][n] = 0;
		}
		
		extMatrix[i][n] = 1;

		ordinary = gaussmethod(extMatrix, n);
		if (ordinary == nullptr)
		{
			cout << "�� ���������� �������� �������." << endl;
			return nullptr;
		}
			for (int j = 0; j < n; j++)
			{
				result[j][i] = ordinary[j];
			}
			
	}

	return result;
}

mytype normVectorUnit(mytype* Vector, int n)
{
	mytype max = 0;
	for (int i = 0; i < n; i++)
	{
		max += abs(Vector[i]);
	}
	return max;
}

mytype normVectorInfinity(mytype* Vector, int n)
{
	mytype max = 0;
	for (int i = 0; i < n; i++)
	{
		if (abs(Vector[i]) > max) max = abs(Vector[i]);
	}
	return max;
}

mytype normMatrixUnit(mytype** Matrix, int n)//���� �������
{
	mytype max = 0;
	mytype *Vector;
	Vector = new mytype[n];
	for (int i = 0; i < n; i++)
	{
		Vector[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Vector[j] += Matrix[i][j];
		}
	}
	max = normVectorInfinity(Vector, n);
	return max;
}

mytype normMatrixInfinity(mytype** Matrix, int n)//���� ������
{
	mytype max = 0;
	mytype *Vector;
	Vector = new mytype[n];
	for (int i = 0; i < n; i++)
	{
		Vector[i] = 0;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Vector[i] += Matrix[i][j];
		}
	}
	max = normVectorInfinity(Vector, n);
	return max;
}

mytype condMatrixUnit(mytype** extMatrix, int n)
{
	mytype cond = 0;
	mytype **Matrix;
	Matrix = new mytype*[n];
	mytype **newreverseMatrix;
	newreverseMatrix = new mytype*[n];

	for (int ii = 0; ii < n; ii++)
	{
		Matrix[ii] = new mytype[n];
		newreverseMatrix[ii] = new mytype[n];

	}

	for (int i = 0; i < n; i++)//��������� �������� �� ������� ������� � Matrix
	{
		for (int j = 0; j < n; j++)
		{
			Matrix[i][j] = extMatrix[i][j];
		}
	}

		newreverseMatrix = reverseMatrix(Matrix, n);
		mytype a = normMatrixUnit(Matrix, n);
		mytype b = normMatrixUnit(newreverseMatrix, n);
		cond = a*b;

		return cond;
}

mytype condMatrixInfinity(mytype** extMatrix, int n)
{
	mytype cond = 0;
	mytype **Matrix;
	Matrix = new mytype*[n];
	mytype **newreverseMatrix;
	newreverseMatrix = new mytype*[n];

	for (int ii = 0; ii < n; ii++)
	{
		Matrix[ii] = new mytype[n];
		newreverseMatrix[ii] = new mytype[n];

	}

	for (int i = 0; i < n; i++)//��������� �������� �� ������� ������� � Matrix
	{
		for (int j = 0; j < n; j++)
		{
			Matrix[i][j] = extMatrix[i][j];
		}
	}

		newreverseMatrix = reverseMatrix(Matrix, n);
		mytype a = normMatrixInfinity(Matrix, n);
		mytype b = normMatrixInfinity(newreverseMatrix, n);
		cond = a*b;

		return cond;
}

mytype estimatecond(mytype** extMatrix, mytype eps, int n)
{
	mytype max = 0;
	mytype *newsolution;
	newsolution = new mytype[n];
	mytype *solution;
	solution = new mytype[n];
	mytype *deltasolution;
	deltasolution = new mytype[n];
	mytype *vectorb;
	vectorb = new mytype[n];

	mytype deltaxrelunit;	
	mytype deltabrelunit;

	mytype deltaxrelinf;
	mytype deltabrelinf;

	mytype normbunit;
	mytype normxunit;

	mytype normbinf;
	mytype normxinf;

    mytype condunit;
	mytype condinf;
	mytype deltaxunit;
	mytype deltaxinf;

	deltasolution = new mytype[n];
	int i;

	for (int j = 0; j < n; j++)//������������ ������ ���������� ������� �� ����������� �������
	{
		vectorb[j] = extMatrix[j][n];
	}

	solution = qrmethod(extMatrix, n, false);

	normbunit = normVectorUnit(vectorb, n);//����� ������� ������ ����� ����� 1
	normbinf = normVectorInfinity(vectorb, n);//����� ������� ������ ����� ����� inf
	deltabrelunit = abs(eps) / normbunit;//������������� ����������� ������� ������ ����� ����� 1
	deltabrelinf = abs(eps) / normbinf;//������������� ����������� ������� ������ ����� ����� inf

		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				extMatrix[j][n] += (1-2*k)*eps;//����������� +0.01 � -0.01
				newsolution = qrmethod(extMatrix, n, false);

				for (int j = 0; j < n; j++)
				{
					deltasolution[j] = newsolution[j] - solution[j];
				}

				//������� ����� ����������� ������� ��� ������������ ���������� ������������� �����������
				deltaxunit = normVectorUnit(deltasolution, n);//����� ������� ������������� ����������� ������� ����� 1
				deltaxinf = normVectorInfinity(deltasolution, n);//����� ������� ������������� ����������� ������� ����� inf
				normxunit = normVectorUnit(solution, n);//����� ������� ����������� ������� ����� 1
				normxinf = normVectorInfinity(solution, n);//����� ������� ����������� ������� ����� inf

				if (abs(normxunit) > 1e-15)
				{
					deltaxrelunit = deltaxunit / normxunit;//������������� ����������� ������� �������
				}

				if (abs(normxinf) > 1e-15)
				{
					deltaxrelinf = deltaxinf / normxinf;//������������� ����������� ������� �������
				}
				else { condinf = deltaxinf *10e+15; }

				condunit = deltaxrelunit / deltabrelunit;
				condinf = deltaxrelinf / deltabrelinf;

				extMatrix[j][n] -= (1 - 2 * k)*eps;
				cout << "������ ����� ��������������� ������ 1 ��� ��������� " << j+1 << "-�� ���������� ������� ������ ����� ��� ��������� �� �� "<< (1 - 2 * k)*eps <<": " << condunit << endl;
				cout << "������ ����� ��������������� ������ inf ��� ��������� " << j+1 << "-�� ���������� ������� ������ ����� ��� ��������� �� �� " << (1 - 2 * k)*eps << ": " << condinf << endl;
				cout << endl;
				if (max < condunit) { max = condunit; }
			}
		}
		cout << "��������� ������ ����� ���������������: ";
	return max;
}

int main()
{
	setlocale(LC_ALL, "rus");
	mytype eps = 0.01;
	mytype cond;
    int k;
	
	mytype **extMatrix = readfromfile("test1.txt", k);
	//double **extMatrix = readfromfile("test2.txt", k);
	//cout << "���� ������� �� �����: " << endl;
	cout << "������� ���� ������� ������: " << endl;
	cout << endl;
	printMatrix(extMatrix, k, k + 1, true);
	cout << endl;


	//cout << "�����: " << endl;
	//double **sqMatrix = getMatrixMinor(extMatrix, 0, 0, k);
	//sqMatrix = getMatrixMinor(sqMatrix, 0, 0, k-1);//2 ������ ������� ����
	//printMatrix(sqMatrix, k-1, k-1, false);
	//cout << "������������: " << endl;
	//cout << determinant(sqMatrix, k - 1) << endl;

	//cout << "���� ������� � ������: " << endl;
	//double **extMatrix1 = readfromscreen(k);
	//printMatrix(extMatrix1, k);

	mytype *solution = gaussmethod(extMatrix, k);
	if (solution != nullptr)
		for (int i = 0; i < k; i++)
		{
			cout << solution[i] << endl;
		}
	cout << endl;
	cout << "����� ������� ������� ����� ";
	mytype discrep1 = discrepancy(extMatrix, solution, k);
	cout << discrep1 << endl;
	cout << endl;
	cout << "����� ��������������� ������� � ������ 1: ";
	mytype cond1 = condMatrixUnit(extMatrix, k);
	cout << cond1;
	cout << endl;
	cout << "����� ��������������� ������� � ������ inf: ";
	mytype cond2 = condMatrixInfinity(extMatrix, k);
	cout << cond2 << endl;
	cout << endl;

	/*cout << "������� ���� ������� ������ � ��������� ������� �������� ��������: " << endl;
	double **extMatrix2 = readfromfile("test2.txt", k);
	printMatrix(extMatrix2, k, k+1, true);
	double *solution1 = gaussmethod(extMatrix2, k);
	if (solution1 != nullptr)
	for (int i = 0; i < k; i++)
	{
	cout << solution1[i] << endl;
	}
	system("pause");
	matrix_destroyer(extMatrix, k);*/
	/*double **extMatrix1 = readfromscreensquare(k);
	double **extMatrix2 = readfromscreensquare(k);
	double **extMatrix3 = multiplyMatrix(extMatrix1, extMatrix2, k);
	printMatrix(extMatrix3, k, k, false);*/

	
	extMatrix = readfromfile("test5.txt", k);
	cout << "������� ���� ������� QR-����������: " << endl;
	cout << endl;
	printMatrix(extMatrix, k, k + 1, true);
	mytype *solution2 = qrmethod(extMatrix, k, true);
	cout << endl;
	if (solution2 != nullptr)
		for (int i = 0; i < k; i++)
		{
			cout << solution2[i] << endl;
		}
	cout << endl;

	cout << "����� ������� ������� ����� ";
	mytype discrep2 = discrepancy(extMatrix, solution2, k);
	cout << discrep2 << endl;
	cout << endl;
	cout << "����� ��������������� ������� � ������ 1: ";
	mytype cond3 = condMatrixUnit(extMatrix, k);
	cout << cond3;
	cout << endl;
	cout << "����� ��������������� ������� � ������ inf: ";
	mytype cond4 = condMatrixInfinity(extMatrix, k);
	cout << cond4 << endl;
	cout << endl;

	cout << "������ ����� ��������������� ������� ���������� ������� � ���������� �������������: ";
	cout << endl;
	cout << endl;
	cond = estimatecond(extMatrix, eps, k);
	cout << cond;
	cout << endl;
	cout << endl;
	cout << "�������� ���������� �������� ������� � ������� ������������ ������� �� �������� � ���: ";
	cout << endl;
	cout << endl;
	mytype **Matrixtest;//��� ��������
	Matrixtest = new mytype*[k];
	mytype **Matrixtestnew;//��� ��������
	Matrixtestnew = new mytype*[k];
	mytype **Matrixtestnewone;//��� ��������
	Matrixtestnewone = new mytype*[k];

		for (int ii = 0; ii < k; ii++)
	{
			Matrixtest[ii] = new mytype[k];
			Matrixtestnew[ii] = new mytype[k];
			Matrixtestnewone[ii] = new mytype[k];
		}
		for (int i = 0; i < k; i++)//��������� �������� �� ������� ������� � ������� r
		{
			for (int j = 0; j < k; j++)
			{
				Matrixtest[i][j] = extMatrix[i][j];
			}
		}
		cout << "������� �: " << endl;
		cout << endl;
		printMatrix(Matrixtest, k, k, false);
		Matrixtestnew = reverseMatrix(Matrixtest, k);
		cout << endl;
		cout << "�������� ������� � �: " << endl;
		cout << endl;
		printMatrix(Matrixtestnew, k, k, false);
		Matrixtestnewone = multiplyMatrix(Matrixtest, Matrixtestnew, k);
		cout << endl;
		cout << "��������� ������������ � � �������� � �: " << endl;
		cout << endl;
		printMatrix(Matrixtestnewone, k, k, false);
		cout << endl;
	
	system("pause");
	return 0;
}
