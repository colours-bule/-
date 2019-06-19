#include <iostream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <tchar.h>
#pragma warning(disable:4996)
using namespace std;

typedef struct     // �������ṹ
{
	double index;
	double x1;
	double y1;
	double x2;
	double y2;
}ImagePoint;
typedef struct      //��õ���ռ丨������XYZ�ṹ
{
	double index;
	double X1;
	double Y1;
	double Z1;
	double X2;
	double Y2;
	double Z2;
}ImageSpacePoint;
typedef struct      //����Ŀ��Ƶ�����ṹ,�ļ������Ƶ㡱Ҫ��
{
	double index;
	double X;
	double Y;
	double Z;
}ControlPoint;
typedef struct
{
	double  index;
	double x;
	double y;
	double z;
	int a;    //���������ʱ���Ƶ�ĵ��
}Scale;
typedef struct //����ģ�͵��������
{
	double index;
	double X1;
	double Y1;
	double Z1;
}Data1;
typedef struct
{
	double index;
	double x1;
	double y1;
	double z1;
}Data2;

void InverseMatrix(double a[], int n)//�������溯�� (������󣬷���������
{
	int *is, *js, i, j, k, l, u, v;
	double d, p;
	is = (int*)malloc(n * sizeof(int));
	js = (int*)malloc(n * sizeof(int));
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				l = i * n + j;
				p = fabs(a[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		if (d + 1.0 == 1.0)
		{
			free(is);
			free(js);
			cout << "�þ��󲻿�����\n";
			return;
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		l = k * n + k;
		a[l] = 1.0 / a[l];
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k * n + j;
				a[u] = a[u] * a[l];
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = i * n + j;
						a[u] = a[u] - a[i * n + k] * a[k * n + j];
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i * n + k;
				a[u] = -a[u] * a[l];
			}
	}  	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i * n + k;
				v = i * n + is[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}
	}
	free(is);
	free(js);
}

int _tmain(int argc, _TCHAR* argv[])//�����в�����һ��������ֵ��ϰ����дΪargc����ʾ��������ʱ�������в����ĸ����ڶ���������ֵ��ϰ����дΪ*argv[]����ʾָ���ַ��������ָ�룬ÿ���ַ�����Ӧһ������
{	//�����ĳ�ֵ	
	double f = 0.153033;	//��Ӱ�����࣬��mΪ��λ
	double bx = 0.2;	   //��Զ����в���bx����mΪ��λ�����Ժ���
	ImagePoint data1[15];
	ImagePoint data2[8];
	ImagePoint data3[11];
	ControlPoint data4[4];
	Data2 data5[5];
	Data1 mddata1[15];
	Data1 mddata2[8];
	Data1 mddata3[11];
	//*************************************************��ȡ�ļ�	
	{
		FILE *fp1;
		fp1 = fopen("C:\\Users\\wangao\\Desktop\\1.txt", "r");    //data1��ʾ·����r��ʾ���ļ���ʽ(Ϊֻ����ʽ������δ����ǽ��ļ�ָ�븳��fp1
		if (fp1 == 0)
		{
			cout << "�޷���ָ���ļ�!\n";
			return 0;
		}
		for (int i = 0; !feof(fp1); i++)   //feof����ļ����ϵ��ļ�������
		{
			fscanf_s(fp1, "%lf %lf %lf %lf %lf ", &data1[i].index, &data1[i].x1, &data1[i].y1, &data1[i].x2, &data1[i].y2);  	//��fp1��ȡ����
			data1[i].x1 = data1[i].x1 / 1000000;		//��΢�ױ�Ϊ��
			data1[i].y1 = data1[i].y1 / 1000000;
			data1[i].x2 = data1[i].x2 / 1000000;
			data1[i].y2 = data1[i].y2 / 1000000;
			cout << "��ȡ�ɹ���\n";
		}
		fclose(fp1);//fp1����ָ��data1
	}
	{
		FILE *fp1;
		fp1 = fopen("C:\\Users\\wangao\\Desktop\\2.txt", "r");
		if (fp1 == 0)
		{
			cout << "�޷���ָ���ļ�!\n";
			return 0;
		}
		for (int i = 0; !feof(fp1); i++)
		{

			fscanf_s(fp1, "%lf %lf %lf %lf %lf ", &data2[i].index, &data2[i].x1, &data2[i].y1, &data2[i].x2, &data2[i].y2);
			data2[i].x1 = data2[i].x1 / 1000000;
			data2[i].y1 = data2[i].y1 / 1000000;
			data2[i].x2 = data2[i].x2 / 1000000;
			data2[i].y2 = data2[i].y2 / 1000000;
			cout << "��ȡ�ɹ���\n";
		}
		fclose(fp1);
	}
	{
		FILE *fp1;
		fp1 = fopen("C:\\Users\\wangao\\Desktop\\3.txt", "r");
		if (fp1 == 0)
		{
			cout << "�޷���ָ���ļ�!\n";
			return 0;
		}
		for (int i = 0; !feof(fp1); i++)
		{
			fscanf_s(fp1, "%lf %lf %lf %lf %lf ", &data3[i].index, &data3[i].x1, &data3[i].y1, &data3[i].x2, &data3[i].y2);
			data3[i].x1 = data3[i].x1 / 1000000;
			data3[i].y1 = data3[i].y1 / 1000000;
			data3[i].x2 = data3[i].x2 / 1000000;
			data3[i].y2 = data3[i].y2 / 1000000;
			cout << "��ȡ�ɹ���\n";
		}
		fclose(fp1);
	}
	{
		FILE *fp1;
		fp1 = fopen("C:\\Users\\wangao\\Desktop\\4.txt", "r");
		if (fp1 == 0)
		{
			cout << "�޷���ָ���ļ�!\n";
			return 0;
		}
		for (int i = 0; !feof(fp1) && i < 4; i++)
		{
			fscanf_s(fp1, "%lf %lf %lf %lf ", &data4[i].index, &data4[i].X, &data4[i].Y, &data4[i].Z);
			cout << "��ȡ�ɹ���\n";
		}
		fclose(fp1);
	}
	{
		FILE *fp1;
		fp1 = fopen("C:\\Users\\wangao\\Desktop\\5.txt", "r");
		if (fp1 == 0)
		{

			cout << "�޷���ָ���ļ�!\n";
			return 0;
		}
		for (int i = 0; !feof(fp1); i++)
		{
			fscanf_s(fp1, "%lf %lf %lf %lf ", &data5[i].index, &data5[i].x1, &data5[i].y1, &data5[i].z1);
			cout << "��ȡ�ɹ���\n";
		}
		fclose(fp1);
	}


	//**********************************************��ȡ�ļ�����  
	//********************************************��Զ���	
	double result1[5] = { 0 };	//��һ����������Զ���Ԫ�صĸ�����
	double result2[5] = { 0 };	//�ڶ�����������Զ���Ԫ�صĸ�����
	double result3[5] = { 0 };	//��������������Զ���Ԫ�صĸ�����
	int i1, i2, i3;
	double U1 = 0, V1 = 0, ��1 = 0, ��1 = 0, ��1 = 0;	 //��ʼ�������Զ���Ԫ��u��v,fai,womiga,��
	double U2 = 0, V2 = 0, ��2 = 0, ��2 = 0, ��2 = 0;
	double U3 = 0, V3 = 0, ��3 = 0, ��3 = 0, ��3 = 0;
	double N11[15], N21[15];	 //15����Ϊ��һ�����������15��ͬ����
	double N12[8], N22[8];		     //8����Ϊ�ڶ������������8��ͬ����
	double N13[11], N23[11];	     //11����Ϊ���������������11��ͬ����
	double a4, a5, a6, b4, b5, b6, c4, c5, c6;
	ImageSpacePoint data12[15];
	double data13[15][5];//����A
	double data14[5][5];
	double R1[5][15];
	double A1[5][15];
	double L1[15];
	double by1;
	double bz1;
	for (i1 = 1;; i1++)
	{
		//R1=E  R2��a1,a2,a3;b1,b2,b3;c1,c2,c3����  		
		//�������Ԫ�س�ʼֵ		  
		U1 = U1 + result1[0];
		V1 = V1 + result1[1];
		��1 = ��1 + result1[2];
		��1 = ��1 + result1[3];
		��1 = ��1 + result1[4];
		by1 = bx * U1;
		bz1 = bx * V1;
		//������ռ�ֱ������תΪ��ռ丨������ϵ�Ĳ���
		double a1, a2, a3, b1, b2, b3, c1, c2, c3;
		a1 = cos(��1)*cos(��1) - sin(��1)*sin(��1)*sin(��1);
		a2 = -cos(��1)*sin(��1) - sin(��1)*sin(��1)*cos(��1);
		a3 = -sin(��1)*cos(��1);
		b1 = cos(��1)*sin(��1);
		b2 = cos(��1)*cos(��1);
		b3 = -sin(��1);
		c1 = sin(��1)*cos(��1) + cos(��1)*sin(��1)*sin(��1);
		c2 = -sin(��1)*sin(��1) + cos(��1)*sin(��1)*cos(��1);
		c3 = cos(��1)*cos(��1);
		//����XYZ��ռ�����ϵ		
		for (int i = 0; i < 15; i++)	//15��ͬ���㣬ѭ��15��	
		{
			//��������������R2��˺����ռ丨������XYZ		
			data12[i].index = data1[i].index;
			data12[i].X1 = data1[i].x1;
			data12[i].Y1 = data1[i].y1;
			data12[i].Z1 = -f;
			data12[i].X2 = a1 * data1[i].x2 + a2 * data1[i].y2 + a3 * (-f);
			data12[i].Y2 = b1 * data1[i].x2 + b2 * data1[i].y2 + b3 * (-f);
			data12[i].Z2 = c1 * data1[i].x2 + c2 * data1[i].y2 + c3 * (-f);
			//����ͶӰϵ��			
			N11[i] = (bx*data12[i].Z2 - bz1 * data12[i].X2) / (data12[i].X1*data12[i].Z2 - data12[i].X2*data12[i].Z1);
			N21[i] = (bx*data12[i].Z1 - bz1 * data12[i].X1) / (data12[i].X1*data12[i].Z2 - data12[i].X2*data12[i].Z1);
			//�������A			
			data13[i][0] = bx;
			data13[i][1] = -(data12[i].Y2 / data12[i].Z2)*bx;
			data13[i][2] = -(data12[i].X2*data12[i].Y2*N21[i]) / data12[i].Z2;
			data13[i][3] = -(data12[i].Z2 + (data12[i].Y2*data12[i].Y2) / data12[i].Z2)*N21[i];
			data13[i][4] = data12[i].X2*N21[i];
			//�������L			
			L1[i] = N11[i] * data12[i].Y1 - N21[i] * data12[i].Y2 - by1;
		}
		//����ת�þ���  	
		for (int i = 0; i < 15; i++)		//����5��ͬ���������	
		{
			for (int j = 0; j < 5; j++)		//������
			{
				R1[j][i] = data13[i][j];
			}
		}
		double t;
		for (int i = 0; i < 5; i++)	//��ATA
		{
			for (int m = 0; m < 5; m++)
			{
				t = 0;
				for (int j = 0; j < 15; j++)
				{
					t = R1[i][j] * data13[j][m] + t;
				}
				data14[i][m] = t;
			}
		}
		InverseMatrix(*data14, 5);  	//��ATA��ת��
		double tt;
		for (int i = 0; i < 5; i++)
		{
			for (int m = 0; m < 15; m++)
			{
				tt = 0;
				for (int j = 0; j < 5; j++)
				{
					tt = data14[i][j] * R1[j][m] + tt;	//ATA��ת����*AT
				}
				A1[i][m] = tt;
			}
		}
		for (int i = 0; i < 5; i++)
		{
			double o = 0;
			for (int j = 0; j < 15; j++)
			{
				o = A1[i][j] * L1[j] + o;	    //ATA��ת��*AT��*L
			}			result1[i] = o;
		}
		if (fabs(result1[0]) < 0.00003 && fabs(result1[1]) < 0.00003 && fabs(result1[2]) < 0.00003 && fabs(result1[3]) < 0.00003 && fabs(result1[4]) < 0.00003)	//��Զ���Ԫ�صĸ�����С��0.00003ʱ��ֹͣѭ����fabs�����ֵ
		{
			��1 = ��1 + result1[2];
			��1 = ��1 + result1[3];
			��1 = ��1 + result1[4];
			a1 = cos(��1)*cos(��1) - sin(��1)*sin(��1)*sin(��1);
			a2 = -cos(��1)*sin(��1) - sin(��1)*sin(��1)*cos(��1);
			a3 = -sin(��1)*cos(��1);
			b1 = cos(��1)*sin(��1);
			b2 = cos(��1)*cos(��1);
			b3 = -sin(��1);
			c1 = sin(��1)*cos(��1) + cos(��1)*sin(��1)*sin(��1);
			c2 = -sin(��1)*sin(��1) + cos(��1)*sin(��1)*cos(��1);
			c3 = cos(��1)*cos(��1);
			a4 = a1;
			a5 = a2;
			a6 = a3;
			b4 = b1;
			b5 = b2;
			b6 = b3;
			c4 = c1;
			c5 = c2;
			c6 = c3;
			break;
		}
	}
	/****************************************///����ڶ�����������ļ�	
	ImageSpacePoint data22[8];
	double data23[8][5];
	double data24[5][5];
	double R2[5][8];
	double A2[5][8];
	double L2[8];
	double by2;
	double bz2;
	for (i2 = 1;; i2++)
	{
		U2 = U2 + result2[0];
		V2 = V2 + result2[1];
		��2 = ��2 + result2[2];
		��2 = ��2 + result2[3];
		��2 = ��2 + result2[4];
		by2 = bx * U2;
		bz2 = bx * V2;

		double  a1, a2, a3, b1, b2, b3, c1, c2, c3;
		a1 = cos(��2)*cos(��2) - sin(��2)*sin(��2)*sin(��2);
		a2 = -cos(��2)*sin(��2) - sin(��2)*sin(��2)*cos(��2);
		a3 = -sin(��2)*cos(��2);
		b1 = cos(��2)*sin(��2);
		b2 = cos(��2)*cos(��2);
		b3 = -sin(��2);
		c1 = sin(��2)*cos(��2) + cos(��2)*sin(��2)*sin(��2);
		c2 = -sin(��2)*sin(��2) + cos(��2)*sin(��2)*cos(��2);
		c3 = cos(��2)*cos(��2);

		for (int i = 0; i < 8; i++)
		{

			data22[i].index = data2[i].index;
			data22[i].X1 = a4 * data2[i].x1 + a5 * data2[i].y1 + a6 * (-f);
			data22[i].X2 = a1 * data2[i].x2 + a2 * data2[i].y2 + a3 * (-f);
			data22[i].Y1 = b4 * data2[i].x1 + b5 * data2[i].y1 + b6 * (-f);
			data22[i].Y2 = b1 * data2[i].x2 + b2 * data2[i].y2 + b3 * (-f);
			data22[i].Z1 = c4 * data2[i].x1 + c5 * data2[i].y1 + c6 * (-f);
			data22[i].Z2 = c1 * data2[i].x2 + c2 * data2[i].y2 + c3 * (-f);

			N12[i] = (bx*data22[i].Z2 - bz2 * data22[i].X2) / (data22[i].X1*data22[i].Z2 - data22[i].X2*data22[i].Z1);
			N22[i] = (bx*data22[i].Z1 - bz2 * data22[i].X1) / (data22[i].X1*data22[i].Z2 - data22[i].X2*data22[i].Z1);

			data23[i][0] = bx;
			data23[i][1] = -(data22[i].Y2 / data22[i].Z2)*bx;
			data23[i][2] = -(data22[i].X2*data22[i].Y2*N22[i]) / data22[i].Z2;
			data23[i][3] = -(data22[i].Z2 + (data22[i].Y2*data22[i].Y2) / data22[i].Z2)*N22[i];
			data23[i][4] = data22[i].X2*N22[i];

			L2[i] = N12[i] * data22[i].Y1 - N22[i] * data22[i].Y2 - by2;
		}  		for (int i = 0; i < 8; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				R2[j][i] = data23[i][j];
			}
		}		double t;
		for (int i = 0; i < 5; i++)
		{
			for (int m = 0; m < 5; m++)
			{
				t = 0;
				for (int j = 0; j < 8; j++)
				{
					t = R2[i][j] * data23[j][m] + t;
				}
				data24[i][m] = t;
			}
		}
		InverseMatrix(*data24, 5);
		//����(ATA��-1AT	
		double tt;
		for (int i = 0; i < 5; i++)
		{
			for (int m = 0; m < 8; m++)
			{
				tt = 0;
				for (int j = 0; j < 5; j++)
				{
					tt = data24[i][j] * R2[j][m] + tt;
				}
				A2[i][m] = tt;
			}
		}
		for (int i = 0; i < 5; i++)
		{
			double o = 0;
			for (int j = 0; j < 8; j++)
			{
				o = A2[i][j] * L2[j] + o;
			}
			result2[i] = o;
		}
		if (fabs(result2[0]) < 0.00003  && fabs(result2[1]) < 0.00003  && fabs(result2[2]) < 0.00003 && fabs(result2[3]) < 0.00003 && fabs(result2[4]) < 0.00003)
		{
			��2 = ��2 + result2[2];
			��2 = ��2 + result2[3];
			��2 = ��2 + result2[4];
			a1 = cos(��2)*cos(��2) - sin(��2)*sin(��2)*sin(��2);
			a2 = -cos(��2)*sin(��2) - sin(��2)*sin(��2)*cos(��2);
			a3 = -sin(��2)*cos(��2);
			b1 = cos(��2)*sin(��2);
			b2 = cos(��2)*cos(��2);
			b3 = -sin(��2);
			c1 = sin(��2)*cos(��2) + cos(��2)*sin(��2)*sin(��2);
			c2 = -sin(��2)*sin(��2) + cos(��2)*sin(��2)*cos(��2);
			c3 = cos(��2)*cos(��2);
			a4 = a1;
			a5 = a2;
			a6 = a3;
			b4 = b1;
			b5 = b2;
			b6 = b3;
			c4 = c1;
			c5 = c2;
			c6 = c3;
			break;
		}
	}
	/****************************************///�����������������ļ�	
	ImageSpacePoint data32[11];
	double data33[11][5];
	double data34[5][5];
	double R3[5][11];
	double A3[5][11];
	double L3[11];
	double by3;
	double bz3;
	for (i3 = 1;; i3++)
	{

		U3 = U3 + result3[0];
		V3 = V3 + result3[1];
		��3 = ��3 + result3[2];
		��3 = ��3 + result3[3];
		��3 = ��3 + result3[4];
		by3 = bx * U3;
		bz3 = bx * V3;

		double a1, a2, a3, b1, b2, b3, c1, c2, c3;
		a1 = cos(��3)*cos(��3) - sin(��3)*sin(��3)*sin(��3);
		a2 = -cos(��3)*sin(��3) - sin(��3)*sin(��3)*cos(��3);
		a3 = -sin(��3)*cos(��3);
		b1 = cos(��3)*sin(��3);
		b2 = cos(��3)*cos(��3);
		b3 = -sin(��3);
		c1 = sin(��3)*cos(��3) + cos(��3)*sin(��3)*sin(��3);
		c2 = -sin(��3)*sin(��3) + cos(��3)*sin(��3)*cos(��3);
		c3 = cos(��3)*cos(��3);

		for (int i = 0; i < 11; i++)
		{

			data32[i].index = data3[i].index;
			data32[i].X1 = a4 * data3[i].x1 + a5 * data3[i].y1 + a6 * (-f);
			data32[i].X2 = a1 * data3[i].x2 + a2 * data3[i].y2 + a3 * (-f);
			data32[i].Y1 = b4 * data3[i].x1 + b5 * data3[i].y1 + b6 * (-f);
			data32[i].Y2 = b1 * data3[i].x2 + b2 * data3[i].y2 + b3 * (-f);
			data32[i].Z1 = c4 * data3[i].x1 + c5 * data3[i].y1 + c6 * (-f);
			data32[i].Z2 = c1 * data3[i].x2 + c2 * data3[i].y2 + c3 * (-f);
			//�����ͶӰϵ��	
			N13[i] = (bx*data32[i].Z2 - bz3 * data32[i].X2) / (data32[i].X1*data32[i].Z2 - data32[i].X2*data32[i].Z1);
			N23[i] = (bx*data32[i].Z1 - bz3 * data32[i].X1) / (data32[i].X1*data32[i].Z2 - data32[i].X2*data32[i].Z1);

			data33[i][0] = bx;
			data33[i][1] = -(data32[i].Y2 / data32[i].Z2)*bx;
			data33[i][2] = -(data32[i].X2*data32[i].Y2*N23[i]) / data32[i].Z2;
			data33[i][3] = -(data32[i].Z2 + (data32[i].Y2*data32[i].Y2) / data32[i].Z2)*N23[i];
			data33[i][4] = data32[i].X2*N23[i];

			L3[i] = N13[i] * data32[i].Y1 - N23[i] * data32[i].Y2 - by3;
		}
		for (int i = 0; i < 11; i++)
		{
			for (int j = 0; j < 5; j++)
			{
				R3[j][i] = data33[i][j];
			}
		}
		double t;
		for (int i = 0; i < 5; i++)
		{
			for (int m = 0; m < 5; m++)
			{
				t = 0;
				for (int j = 0; j < 11; j++)
				{
					t = R3[i][j] * data33[j][m] + t;
				}			data34[i][m] = t;
			}
		}
		InverseMatrix(*data34, 5);
		//����(ATA��-1AT  	
		double tt;
		for (int i = 0; i < 5; i++)
		{
			for (int m = 0; m < 11; m++)
			{
				tt = 0;
				for (int j = 0; j < 5; j++)
				{
					tt = data34[i][j] * R3[j][m] + tt;
				}			A3[i][m] = tt;
			}
		}
		for (int i = 0; i < 5; i++)
		{
			double o = 0;
			for (int j = 0; j < 11; j++)
			{
				o = A3[i][j] * L3[j] + o;
			}		result3[i] = o;
		}
		if (fabs(result3[0]) < 0.000003 && fabs(result3[1]) < 0.000003  &&  fabs(result3[2]) < 0.000003  && fabs(result3[3]) < 0.000003  && fabs(result3[4]) < 0.000003)
		{
			break;
		}
	}//***********************************************************************��Զ������
	//***********************************************************************��ʼģ�����ӣ�����������֮��ı����߹黯ϵ��k)
	//����ģ�͵��ڵ����������ϵ�е����꣨û�о��������ģ� 
	for (int i = 0; i < 15; i++)//���һ
	{
		mddata1[i].index = data1[i].index;
		mddata1[i].X1 = data12[i].X1*N11[i];
		mddata1[i].Z1 = data12[i].Z1*N11[i];
		mddata1[i].Y1 = (data12[i].Y1*N11[i] + data12[i].Y2*N21[i] + by1) / 2;//Yȡƽ��ֵ�ǿ��ǵ����������Ӳ��Ӱ�죬�������Լ���Ӱ��
	}
	for (int i = 0; i < 8; i++)//��Զ�
	{
		mddata2[i].index = data2[i].index;
		mddata2[i].X1 = data22[i].X1*N12[i];
		mddata2[i].Z1 = data22[i].Z1*N12[i];
		mddata2[i].Y1 = (data22[i].Y1*N12[i] + data22[i].Y2*N22[i] + by2) / 2;
	}
	for (int i = 0; i < 11; i++)  //�����
	{
		mddata3[i].index = data3[i].index;
		mddata3[i].X1 = data32[i].X1*N13[i];
		mddata3[i].Z1 = data32[i].Z1*N13[i];
		mddata3[i].Y1 = (data32[i].Y1*N13[i] + data32[i].Y2*N23[i] + by3) / 2;
	}
	// ����ÿ�����������������ı����߱任ϵ��k
	double K1[3] = { 0 }, k1 = 0;
	int count1 = 0;
	double K2[3] = { 0 }, k2 = 0;
	int count2 = 0;
	for (int i = 0; i < 15; i++)//��ģ��һ��ģ�Ͷ�Ϊ��׼
	{
		for (int j = 0; j < 8; j++)
		{
			if (data1[i].index == data2[j].index)
			{
				K2[count2] = (mddata1[i].Z1 - bz1) / (mddata2[j].Z1);
				k1 = k1 + K2[count2];
				count2++;
			}
		}
	}
	for (int i = 0; i < 8; i++)//��ģ�Ͷ���ģ����Ϊ��׼
	{
		for (int j = 0; j < 11; j++)
		{
			if (data22[i].index == data32[j].index)
			{
				K1[count1] = (mddata2[i].Z1 - bz2) / (mddata3[j].Z1);
				k2 = k2 + K1[count1];
				count1++;
			}
		}
	}
	k1 = k1 / 3;
	k2 = k2 / 3;
	//******************************************************ģ�����ӽ���
	//******************************************************������Ӱ������
	//����Ƶ��Ӧ���������
	Scale P1[10] = { 0 };
	int point1 = 0;
	for (int i = 0; i < 15; i++)//�����1��Ѱ�Ҷ�Ӧ�Ŀ��Ƶ�
	{
		for (int j = 0; j < 4; j++)
		{
			if (mddata1[i].index == data4[j].index)
			{
				P1[point1].index = mddata1[i].index;
				P1[point1].x = mddata1[i].X1;
				P1[point1].y = mddata1[i].Y1;
				P1[point1].z = mddata1[i].Z1;
				point1++;
			}
		}
	}
	Scale P2[10] = { 0 };
	int point2 = 0;
	for (int i = 0; i < 11; i++)//�����3��Ѱ�Ҷ�Ӧ�Ŀ��Ƶ㣨�������2�е�ģ�͵�ֱ������1��3�ж��У��������������2���ٴν��в��ң�
	{
		for (int j = 0; j < 4; j++)
		{
			if (mddata3[i].index == data4[j].index)
			{
				P2[point2].index = mddata3[i].index;
				P2[point2].x = mddata3[i].X1;
				P2[point2].y = mddata3[i].Y1;
				P2[point2].z = mddata3[i].Z1;
				point2++;
			}
		}
	}
	double blc[3]; //������
	double length[4]; //���Ƶ�������� 
	//����Ƶ��Ӧ������ͼ�Ͼ���
	int i = 1;
	blc[0] = sqrt((P1[i].x - P1[i + 1].x)*(P1[i].x - P1[i + 1].x) + (P1[i].y - P1[i + 1].y)*(P1[i].y - P1[i + 1].y));
	blc[1] = sqrt((P1[i].x - P1[i - 1].x)*(P1[i].x - P1[i - 1].x) + (P1[i].y - P1[i - 1].y)*(P1[i].y - P1[i - 1].y));
	blc[2] = sqrt((P1[i - 1].x - P1[i + 1].x)*(P1[i - 1].x - P1[i + 1].x) + (P1[i - 1].y - P1[i + 1].y)*(P1[i - 1].y - P1[i + 1].y));
	//����Ƶ��ʵ�ʾ���
	length[0] = sqrt((data4[0].X - data4[1].X)*(data4[0].X - data4[1].X) + (data4[0].Y - data4[1].Y)*(data4[0].Y - data4[1].Y));
	length[1] = sqrt((data4[2].X - data4[1].X)*(data4[2].X - data4[1].X) + (data4[2].Y - data4[1].Y)*(data4[2].Y - data4[1].Y));
	length[2] = sqrt((data4[3].X - data4[1].X)*(data4[3].X - data4[1].X) + (data4[3].Y - data4[1].Y)*(data4[3].Y - data4[1].Y));
	length[3] = sqrt((data4[0].X - data4[2].X)*(data4[0].X - data4[2].X) + (data4[0].Y - data4[2].Y)*(data4[0].Y - data4[2].Y));
	//��ʼ���������
	double  m[3], m0;
	m[0] = length[1] / blc[0];
	m[1] = length[0] / blc[1];
	m[2] = length[3] / blc[2];
	m0 = (m[0] + m[1] + m[2]) / 3;
	//********************************************************������Ӱ�����߽���  

	//********************************************************�������ɺ�����
	//�����ģ����վ��Ӱ��������
	double Xps[3], Yps[3], Zps[3];//��վ����Ӱ��������(Xps,Yps,Zps)
	for (int i = 0; i < 3; i++)
	{
		if (i == 0)	//��һ��ģ������վ
		{
			Xps[0] = 0;
			Yps[0] = 0;
			Zps[0] = m0 * f;
		}
		else if (i == 1)//��һ��ģ������վ
		{
			Xps[i] = Xps[i - 1] + m0 * bx;
			Yps[i] = Yps[i - 1] + m0 * by1;
			Zps[i] = Zps[i - 1] + m0 * bz1;
		}
		else        //��һ��ģ���Ժ����ģ�͵�����վ���������
		{
			Xps[i] = Xps[i - 1] + m0 * k1* bx;
			Yps[i] = Yps[i - 1] + m0 * k1*by2;
			Zps[i] = Zps[i - 1] + m0 * k1*bz2;
		}
	}

	//�������е�ģ�͵��������
	for (int i = 0; i < 15; i++)//���1
	{
		mddata1[i].X1 = mddata1[i].X1*m0;
		mddata1[i].Z1 = mddata1[i].Z1*m0 + m0 * f;
		mddata1[i].Y1 = mddata1[i].Y1*m0;
	}
	for (int i = 0; i < 8; i++)//���2
	{
		mddata2[i].X1 = Xps[1] + k1 * mddata2[i].X1*m0;
		mddata2[i].Y1 = (Yps[1] + k1 * m0*N12[i] * data22[i].Y1 + k1 * m0*N22[i] * data22[i].Y2 + Yps[2]) / 2;
		mddata2[i].Z1 = Zps[1] + k1 * mddata2[i].Z1*m0;
	}
	for (int i = 0; i < 11; i++)//���3
	{
		mddata3[i].X1 = Xps[2] + k2 * k1*mddata3[i].X1*m0;
		mddata3[i].Y1 = (Yps[2] + k2 * k1*m0*N13[i] * data32[i].Y1 + k2 * k1*m0*N23[i] * data32[i].Y2 + Yps[2] + k2 * k1*m0*by3) / 2;
		mddata3[i].Z1 = Zps[2] + k2 * k1*mddata3[i].Z1*m0;
	}
	//���е�ģ�͵�����������
	//*********************************************************�������������

	//*********************************************************��ʼ����ģ�͵ľ��Զ���
	// ���������ֵ�����õ�һ�����Ƶ������һ�����Ƶ���ɣ�
	double dXt, dYt, dXp, dYp;
	dXt = data4[0].X - data4[3].X;//�����������
	dYt = data4[0].Y - data4[3].Y;
	dXp = mddata1[0].X1 - mddata3[8].X1;//�����������
	dYp = mddata1[0].Y1 - mddata3[8].Y1;
	//������Զ������������a b T
	double a, b, T;
	a = (dXp*dYt + dYp * dXt) / ((dXt*dXt) + (dYt*dYt));
	b = (dXp*dXt - dYp * dYt) / ((dXt*dXt) + (dYt*dYt));
	T = sqrt(a*a + b * b);


	//������Ƶ�任�ɵĵ�����Ӱ�������� 
	ControlPoint data43[4];//���Ƶ�ĵ�����Ӱ��������
	for (int i = 0; i < 4; i++)
	{
		data43[i].X = b * (data4[i].X - data4[0].X) + (data4[i].Y - data4[0].Y)*a;
		data43[i].Y = a * (data4[i].X - data4[0].X) + (data4[i].Y - data4[0].Y)*(-b);
		data43[i].Z = T * (data4[i].Z - data4[0].Z);
	}
	//��Ӱ��������Xp
	Scale P3[10] = { 0 };
	point1 = 0;
	for (int i = 0; i < 15; i++)//�����1��Ѱ�Ҷ�Ӧ�Ŀ��Ƶ�
	{
		for (int j = 0; j < 4; j++)
		{
			if (mddata1[i].index == data4[j].index)
			{
				P3[point1].index = mddata1[i].index;
				P3[point1].x = mddata1[i].X1;
				P3[point1].y = mddata1[i].Y1;
				P3[point1].z = mddata1[i].Z1;
				point1++;
			}
		}
	}
	Scale P4[10] = { 0 };
	point2 = 0;
	for (int i = 0; i < 11; i++)//�����3��Ѱ�Ҷ�Ӧ�Ŀ��Ƶ�
	{
		for (int j = 0; j < 4; j++)
		{
			if (mddata3[i].index == data4[j].index)
			{
				P4[point2].index = mddata3[i].index;
				P4[point2].x = mddata3[i].X1;
				P4[point2].y = mddata3[i].Y1;
				P4[point2].z = mddata3[i].Z1;
				point2++;
			}
		}
	}

	//���������Ӱ�����������Ӱ��������
	ControlPoint tp[4], p[4];
	for (int i = 0; i < 4; i++)
	{
		tp[i].X = data43[i].X;
		tp[i].Y = data43[i].Y;
		tp[i].Z = data43[i].Z;
	}
	for (int i = 0; i < 3; i++)
	{
		p[i].X = P3[i].x;
		p[i].Y = P3[i].y;
		p[i].Z = P3[i].z;
	}
	p[3].X = P4[0].x;
	p[3].Y = P4[0].y;
	p[3].Z = P4[0].z;


	//������Զ���Ԫ��
	double �� = 1, X0 = 0, Y0 = 0, Z0 = 0, �� = 0, �� = 0, �� = 0;
	double result4[7] = { 0 };
	double L4[12] = { 0 };
	double A[12][7] = { 0 };
	double AT[7][12];
	double ATA[7][7];
	double AA[7][12];
	int ii = 0;
	for (ii = 1; ; ii++)
	{	//R1=E  R2��a1,a2,a3;b1,b2,b3;c1,c2,c3����  
		//�����ʼֵ	  		
		X0 = X0 + result4[0];
		Y0 = Y0 + result4[1];
		Z0 = Z0 + result4[2];
		�� = �� + result4[3];
		�� = �� + result4[4];
		�� = �� + result4[5];
		�� = �� + result4[6];
		//������ռ�ֱ������תΪ��ռ丨������ϵ�Ĳ���		
		double  a1, a2, a3, b1, b2, b3, c1, c2, c3;
		a1 = cos(��)*cos(��) - sin(��)*sin(��)*sin(��);
		a2 = -cos(��)*sin(��) - sin(��)*sin(��)*cos(��);
		a3 = -sin(��)*cos(��);
		b1 = cos(��)*sin(��);
		b2 = cos(��)*cos(��);
		b3 = -sin(��);
		c1 = sin(��)*cos(��) + cos(��)*sin(��)*sin(��);
		c2 = -sin(��)*sin(��) + cos(��)*sin(��)*cos(��);
		c3 = cos(��)*cos(��);
		//���㳣�������L	
		for (int i = 0; i < 4; i++)
		{
			int j = i * 3;
			L4[j] = tp[i].X - �� * (a1*p[i].X + a2 * p[i].Y + a3 * p[i].Z) - X0;
			L4[j + 1] = tp[i].Y - �� * (b1*p[i].X + b2 * p[i].Y + b3 * p[i].Z) - Y0;
			L4[j + 2] = tp[i].Z - �� * (c1*p[i].X + c2 * p[i].Y + c3 * p[i].Z) - Z0;
		}
		//�������A
		for (int i = 0; i < 12; i++)
		{
			for (int j = 0; j < 7; j++)
			{
				if (i % 3 == j)
				{
					A[i][j] = 1;
				}
				else if ((i == 0 || i == 3 || i == 6 || i == 9) && j == 3)
				{
					A[i][j] = p[i / 3].X;
				}
				else if ((i == 2 || i == 5 || i == 8 || i == 11) && j == 4)
				{
					A[i][j] = p[i / 3].X;
				}
				else if ((i == 1 || i == 4 || i == 7 || i == 10) && j == 6)
				{
					A[i][j] = p[i / 3].X;
				}
				else if ((i == 1 || i == 4 || i == 7 || i == 10) && j == 3)
				{
					A[i][j] = p[i / 3].Y;
				}
				else if ((i == 2 || i == 5 || i == 8 || i == 11) && j == 5)
				{
					A[i][j] = p[i / 3].Y;
				}
				else if ((i == 0 || i == 3 || i == 6 || i == 9) && j == 6)
				{
					A[i][j] = -p[i / 3].Y;
				}
				else if ((i == 2 || i == 5 || i == 8 || i == 11) && j == 3)
				{
					A[i][j] = p[i / 3].Z;
				}
				else if ((i == 1 || i == 4 || i == 7 || i == 10) && j == 5)
				{
					A[i][j] = -p[i / 3].Z;
				}
				else if ((i == 0 || i == 3 || i == 6 || i == 9) && j == 4)
				{
					A[i][j] = -p[i / 3].Z;
				}
			}
		}
		//�������A��ת�þ���
		for (int i = 0; i < 12; i++)
		{
			for (int j = 0; j < 7; j++)
			{
				AT[j][i] = A[i][j];
			}
		}
		//�������A��ת�ú�A�ĳ˻�	
		double t;
		for (int i = 0; i < 7; i++)
		{
			for (int m = 0; m < 7; m++)
			{
				t = 0;
				for (int j = 0; j < 12; j++)
				{
					t = AT[i][j] * A[j][m] + t;
				}
				ATA[i][m] = t;
			}
		}
		//����	
		InverseMatrix(*ATA, 7);
		//����(ATA��-1AT	
		double tt;
		for (int i = 0; i < 7; i++)
		{
			for (int m = 0; m < 12; m++)
			{
				tt = 0;
				for (int j = 0; j < 7; j++)
				{
					tt = ATA[i][j] * AT[j][m] + tt;
				}
				AA[i][m] = tt;
			}
		}
		//�������ո������
		for (int i = 0; i < 7; i++)
		{
			double o = 0;
			for (int j = 0; j < 12; j++)
			{
				o = AA[i][j] * L4[j] + o;
			}
			result4[i] = o;
		}
		if (fabs(result4[0]) < 0.000001  &&  fabs(result4[1]) < 0.0000001  &&  fabs(result4[2]) < 0.0000001 && fabs(result4[3]) < 0.0000001 && fabs(result4[4]) < 0.0000001  && fabs(result4[5]) < 0.0000001 && fabs(result4[6]) < 0.0000001)
		{
			break;
		}
	}
	//������Զ���Ԫ�ؽ���
	//*****************************************************����ģ�͵ľ��Զ������

	// ����Ӱ��������任Ϊ������Ӱ��������     
	X0 = X0 + result4[0];
	Y0 = Y0 + result4[1];
	Z0 = Z0 + result4[2];
	�� = �� + result4[3];
	�� = �� + result4[4];
	�� = �� + result4[5];
	�� = �� + result4[6];
	ControlPoint Down1[15], Down2[8], Down3[11];//�������������Ӱ�������� 
	for (int i = 0; i < 15; i++)//���1
	{
		Down1[i].index = mddata1[i].index;
		Down1[i].X = �� * (mddata1[i].X1 - �� * mddata1[i].Y1 - �� * mddata1[i].Z1) + X0;
		Down1[i].Y = �� * (��*mddata1[i].X1 + mddata1[i].Y1 - �� * mddata1[i].Z1) + Y0;
		Down1[i].Z = �� * (��*mddata1[i].X1 + �� * mddata1[i].Y1 + mddata1[i].Z1) + Z0;
	}
	for (int i = 0; i < 8; i++)//���2
	{
		Down2[i].index = mddata2[i].index;
		Down2[i].X = �� * (mddata2[i].X1 - �� * mddata2[i].Y1 - �� * mddata2[i].Z1) + X0;
		Down2[i].Y = �� * (��*mddata2[i].X1 + mddata2[i].Y1 - �� * mddata2[i].Z1) + Y0;
		Down2[i].Z = �� * (��*mddata2[i].X1 + �� * mddata2[i].Y1 + mddata2[i].Z1) + Z0;
	}
	for (int i = 0; i < 11; i++)//���3
	{
		Down3[i].index = mddata3[i].index;
		Down3[i].X = �� * (mddata3[i].X1 - �� * mddata3[i].Y1 - �� * mddata3[i].Z1) + X0;
		Down3[i].Y = �� * (��*mddata3[i].X1 + mddata3[i].Y1 - �� * mddata3[i].Z1) + Y0;
		Down3[i].Z = �� * (��*mddata3[i].X1 + �� * mddata3[i].Y1 + mddata3[i].Z1) + Z0;
	}
	//�������յĴ������
	ControlPoint Ground1[15], Ground2[8], Ground3[11];//����������
	for (int i = 0; i < 15; i++)
	{
		Ground1[i].index = Down1[i].index;
		Ground1[i].X = (b*Down1[i].X + a * Down1[i].Y) / T / T + data4[0].X;
		Ground1[i].Y = (a*Down1[i].X - b * Down1[i].Y) / T / T + data4[0].Y;
		Ground1[i].Z = (Down1[i].Z) / T + data4[0].Z;
	}
	for (int i = 0; i < 8; i++)
	{
		Ground2[i].index = Down2[i].index;
		Ground2[i].X = (b*Down2[i].X + a * Down2[i].Y) / T / T + data4[0].X;
		Ground2[i].Y = (a*Down2[i].X - b * Down2[i].Y) / T / T + data4[0].Y;
		Ground2[i].Z = (Down2[i].Z) / T + data4[0].Z;
	}
	for (int i = 0; i < 11; i++)
	{
		Ground3[i].index = Down3[i].index;
		Ground3[i].X = (b*Down3[i].X + a * Down3[i].Y) / T / T + data4[0].X;
		Ground3[i].Y = (a*Down3[i].X - b * Down3[i].Y) / T / T + data4[0].Y;
		Ground3[i].Z = (Down3[i].Z) / T + data4[0].Z;
	}

	//�������������Ĳ�ֵ	  		
	double  D[5][4] = { 0 };
	for (i = 0; i < 5; i++)
	{
		for (int j = 0; j < 15; j++)
		{
			if (data5[i].index == Ground1[j].index)
			{
				D[i][0] = data5[i].index;
				D[i][1] = data5[i].x1 - Ground1[j].X;
				D[i][2] = data5[i].y1 - Ground1[j].Y;
				D[i][3] = data5[i].z1 - Ground1[j].Z;
				break;
			}
		}
		for (int j = 0; j < 8; j++)
		{
			if (data5[i].index == Ground2[j].index)
			{
				D[i][0] = data5[i].index;
				D[i][1] = data5[i].x1 - Ground2[j].X;
				D[i][2] = data5[i].y1 - Ground2[j].Y;
				D[i][3] = data5[i].z1 - Ground2[j].Z;
				break;
			}
		}
		for (int j = 0; j < 11; j++)
		{
			if (data5[i].index == Ground3[j].index)
			{
				D[i][0] = data5[i].index;
				D[i][1] = data5[i].x1 - Ground3[j].X;
				D[i][2] = data5[i].y1 - Ground3[j].Y;
				D[i][3] = data5[i].z1 - Ground3[j].Z;
				break;
			}
		}
	}

	//****************************************************************д���ļ�
	FILE *fp2;
	fp2 = fopen("C:\\Users\\wangao\\Desktop\\Result.txt", "w");
	fprintf(fp2, " ��������ϸ���ģ�͵�Ĵ������ \n");
	fprintf(fp2, " X��m)                   Y(m)                     Z(m) \n");
	fprintf(fp2, " ���һ \n");
	for (int j = 0; j < 15; j++)
	{
		fprintf(fp2, " %f   %f    %f \n", Ground1[j].X, Ground1[j].Y, Ground1[j].Z);
	}
	fprintf(fp2, "\n");
	fprintf(fp2, " ��Զ� \n");
	for (int j = 0; j < 8; j++)
	{
		fprintf(fp2, " %f   %f    %f\n ", Ground2[j].X, Ground2[j].Y, Ground2[j].Z);
	}
	fprintf(fp2, "\n");
	fprintf(fp2, " ����� \n");
	for (int j = 0; j < 11; j++)
	{
		fprintf(fp2, " %f   %f    %f \n", Ground3[j].X, Ground3[j].Y, Ground3[j].Z);
	}
	fprintf(fp2, "\n �����Ĳ�ֵ \n");
	fprintf(fp2, "dX               dY               dZ \n");
	for (int j = 0; j < 5; j++)
	{
		fprintf(fp2, " %f   %f    %f  \n", D[j][1], D[j][2], D[j][3]);
	}
	fprintf(fp2, "\n");
	fclose(fp2);
	//*****************************************************************************д���ļ�����
	return 0;
}

