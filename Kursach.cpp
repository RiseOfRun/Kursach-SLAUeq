﻿#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <iomanip>
#include <fstream>
// M = l * C

using namespace std;
typedef vector<vector<double>> Matrix;

struct MatrixProf
{
	int size;
	vector<double> DI;
	vector<double> AL;
	vector<double> AU;
	vector<int> IA;
	vector<int> JA;
};

struct AuxVectors
{
	vector<double> Ax;
	vector<double> r;
	vector<double> z;
	vector<double> p;
	vector<double> LU;
	vector<double> temp;
};


class Net
{
public:
	Net()
	{

	}

	Net(fstream & nodes, fstream & elements, fstream & fields, fstream & condi1, fstream & condi2, fstream & condi3)
	{
		double x, y;
		int a, b, c, field,type;

		while (nodes >> x >> y)
		{
			Node.push_back({ x,y });
		}

		while (elements >>a >>b >>c)
		{
			Elements.push_back({ a,b,c });
		}

		while (fields >> field)
		{
			this->fields.push_back(field);
		}

		while (condi1 >> a >> type)
		{
			this->firstCondi.push_back({a,type});
		}

		while (condi2 >> a >> b >> type)
		{
			this->SecondCondi.push_back({ a, b,type });
		}

		while (condi3 >> a >> b  >> type >> field)
		{
			this->ThirdCondi.push_back({ a,b,type, field });
		}
	}

	void SaveNet(fstream & nodes, fstream & elements, fstream & fields)
	{
		int length = Node.size();
		for (size_t i = 0; i < length; i++)
		{
			nodes << Node[i][0] << " " << Node[i][1] << "\n";
		}

		length = this->Elements.size();
		for (size_t i = 0; i < length; i++)
		{
			elements << this->Elements[i][0] << " " << this->Elements[i][1]<<" "<< this->Elements[i][2] << "\n";
			fields << this->fields[i] << "\n";
		}
	}
	vector<vector<double>> Node;
	vector<vector<int>> Elements;
	vector<int> fields;

	vector<vector<int>> firstCondi;
	vector<vector<int>> SecondCondi;
	vector<vector<int>> ThirdCondi;
	//добавление информации о подобластях
	void DevideBy2Fields()
	{
		fields = vector<int>(Elements.size());
		int middle = fields.size() / 2;
		for (int i = middle; i < fields.size(); i++)
		{
			fields[i] = 1;
		}
	}
	//
	void AddCondi(int nx, int ny)
	{
		for (int i = 0; i < nx; i++)
		{
			firstCondi.push_back({ i,0 });
		}
	}

	//построение сетки на треугольниках
	void BuildNet(double xmin, double xmax, double ymin, double ymax, int nx, int ny)
	{
		
		double hx = (xmax - xmin) / nx;
		double hy = (ymax - ymin) / ny;
		Node = vector<vector<double>>((nx + 1) * (ny + 1));
		Node[0] = vector<double>{ xmin, ymin };

		for (int i = 0; i < ny; i++)
		{
			double y = ymin + i * hy;
			for (int j = 0; j < nx; j++)
			{
				double x = xmin + j * hx;
				Node[i * (nx + 1) + j + 1] = { x + hx, y };
				Node[(i + 1) * (nx + 1) + j] = { x,y + hy };
				Elements.push_back({ j + i * (nx + 1),j + 1 + i * (nx + 1), j + (nx + 1) * (i + 1) });
			}
		}
		Node[Node.size() - 1] = { xmax,ymax };

		for (int i = ny; i > 0; i--)
		{
			for (int j = nx; j > 0; j--)
			{
				Elements.push_back({ j + i * (nx + 1) - nx - 1,j - 1 + i * (nx + 1), j + i * (nx + 1) });
			}
		}
		int length = Elements.size();
		vector<vector<int>> Elementstmp(length);
		for (int j = 0, i = 0; i < length; j++, i += 2)
		{
			Elementstmp[i] = Elements[j];
		}

		for (int i = 1, j = length - 1; i < length; i += 2, j--)
		{
			Elementstmp[i] = Elements[j];
		}

		Elements = Elementstmp;
		fields.resize(Elements.size());
		//разбиение на подобласти
		//DevideBy2Fields();
	}
private:
};


class Eq
{
public:

	Net TheNet;
	vector<double> b;
	MatrixProf AProf;
	MatrixProf LU;
	Matrix A;
	vector<double>q;
	double BigEl;

	void PrintPlot(Matrix& A)
	{
		int length = A.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				cout << A[i][j] << " \t";
			}
			cout << "\n";
		}
	}

	//параметры
	double U(double x, double y)
	{
		return x+y;
	}

	Eq()
	{
		TheNet.BuildNet(0, 2, 0, 2, 2, 2);
		A = Matrix(TheNet.Node.size());
		b = vector<double>(TheNet.Node.size());
	}

	Eq(Net net)
	{
		TheNet = net;
		A = Matrix(TheNet.Node.size());
		for (int i = 0; i < A.size(); i++)
		{
			A[i] = vector<double>(A.size());
		}
		b = vector<double>(TheNet.Node.size());
		q = vector<double>(TheNet.Node.size());
	}
	
	//построение матрицы G
	vector<vector<double>> BuildG(vector<vector<double>>& D_1, double DetD, vector<int>& el, int field) {
		vector<vector<double>> G(3);
		double multix = abs(DetD) / 2;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{

				G[i].push_back(Lambda(TheNet.Node[0], field) * multix * (D_1[i][1] * D_1[j][1] + D_1[i][2] * D_1[j][2])); // Lambda = const;
			}
		}
		return G;
	}
	//разложение коэф дифузии по линейным базисным функциям
	vector<vector<double>> BuildGDecomposeLinalL(vector<vector<double>>& D_1, double DetD, vector<int>& el, int field) {
		vector<vector<double>> G(3);
		double multix = abs(DetD) / 6.;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double l1 = Lambda(TheNet.Node[el[0]], field);
				double l2 = Lambda(TheNet.Node[el[1]], field);
				double l3 = Lambda(TheNet.Node[el[2]], field);

				double sumL = Lambda(TheNet.Node[el[0]], field) +
					Lambda(TheNet.Node[el[1]], field) +
					Lambda(TheNet.Node[el[2]], field);
				G[i].push_back(sumL * multix * (D_1[i][1] * D_1[j][1] + D_1[i][2] * D_1[j][2])); // Lambda = const;
			}
		}
		return G;
	}
	double findMax(double x1, double x2)
	{
		if (x1 > x2)
			return x1;
		else
			return x2;
	}
	//потроение матрицы С, где M = gamma * C
	Matrix BuildC(double DetD)
	{
		Matrix M = Matrix{ {2,1,1 }, { 1,2,1 }, { 1,1,2 } };
		double mult = abs(DetD) / 24;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				M[i][j] *= mult;
			}
		}

		return M;
	}
	//умножение матрицы на вектор
	vector<double> MVecMult(Matrix& A, vector<double>& b)
	{
		vector<double> result(A.size());
		int length = A.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				result[i] += A[i][j] * b[j];
			}
		}
		return result;
	}

	
	void BuildGlobalProf()
	{
		BuildProfile();

		for (int i = 0; i < TheNet.Elements.size(); i++)
		{
			vector<int> element = TheNet.Elements[i];
			int field = TheNet.fields[i];
			BuildLocal(element, field);
			//PrintPlot(A);
			cout << "\n\n\n";
		}
	}

	Matrix BuildLocal(vector<int>& el, int field)
	{
		double x1 = TheNet.Node[el[0]][0];
		double x2 = TheNet.Node[el[1]][0];
		double x3 = TheNet.Node[el[2]][0];
		double y1 = TheNet.Node[el[0]][1];
		double y2 = TheNet.Node[el[1]][1];
		double y3 = TheNet.Node[el[2]][1];

		vector<vector<double>> D{
		vector<double>{1,1,1},
		vector<double> {x1,x2,x3},
		vector<double> {y1,y2,y3}
		};

		double DetD = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);

		vector<vector<double>> D_1{
		vector<double> {x2* y3 - x3 * y2, y2 - y3, x3 - x2},
		vector<double> {x3* y1 - x1 * y3, y3 - y1, x1 - x3},
		vector<double> {x1* y2 - x2 * y1, y1 - y2, x2 - x1}
		};

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				D_1[i][j] /= DetD;
			}
		}
		vector<vector<double>> G = BuildGDecomposeLinalL(D_1, DetD, el, field); //LAMBDA RAZLOJENIE ToRemove
		Matrix M = BuildC(DetD);
		vector<double> f = { F(x1,y1,field),F(x2,y2,field),F(x3,y3,field) };
		vector<double> b = MVecMult(M, f);
		int length = 3;
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				G[i][j] += M[i][j] * Gamma(field);
			}
		}
		ToGLobalProf(G, b, el);
		ToGlobalPlot(G, b, el);
		//double mesP = 1 / 2 * abs(DetD); //РїР»РѕС‰Р°РґСЊ С‚СЂРµСѓРіРѕР»СЊРЅРѕРіРѕ РєСѓСЃРѕС‡РєР° РѕР±Р»Р°СЃС‚Рё
		//ToGlobalPlot(G, b, el);
		return G;
	}

	void BuildGlobalPlot()
	{
		for (int i = 0; i < TheNet.Elements.size(); i++)
		{
			Matrix test = BuildLocal(TheNet.Elements[i], TheNet.fields[i]);
		}
	}

	void BuildProfile()
	{
		vector<vector<int>> profile(TheNet.Node.size());

		for (int i = 0; i < TheNet.Elements.size(); i++)
		{
			for (int j = 1; j < 3; j++)
			{
				for (int k = 0; k < j; k++)
				{
					int current = TheNet.Elements[i][j];
					int node = TheNet.Elements[i][k];
					if (!count(profile[current].begin(), profile[current].end(), node))
					{
						if (profile[current].size() != 0 && profile[current][profile[current].size() - 1] > node)
						{
							for (int l = 0; l < profile[current].size(); l++)
							{
								if (node < profile[current][l])
								{
									profile[current].insert(profile[current].begin() + l, node);
									break;
								}
							}
						}
						else
						{
							profile[current].push_back(node);
						}
					}
				}
			}
		}

		AProf.IA.push_back(1);
		int count = 0;
		for (int i = 1; i < TheNet.Node.size(); i++)
		{
			AProf.IA.push_back(AProf.IA[i - 1] + count);
			count = 0;
			for (int j = 0; j < profile[i].size(); j++)
			{
				AProf.JA.push_back(profile[i][j]);
				count++;
			}
		}
		AProf.IA.push_back(AProf.IA[AProf.IA.size() - 1] + count);
		AProf.AL = vector<double>(AProf.IA[AProf.IA.size() - 1] - 1);
		AProf.AU = vector<double>(AProf.IA[AProf.IA.size() - 1] - 1);
		AProf.DI = vector<double>(TheNet.Node.size());
		AProf.size = AProf.DI.size();
	}

	void AddThirdCondi()
	{
		int length = TheNet.ThirdCondi.size();
		for (int i = 0; i < length; i++)
		{
			vector<int> Edge = TheNet.ThirdCondi[i];
			double x1 = TheNet.Node[Edge[0]][0];
			double y1 = TheNet.Node[Edge[0]][1];
			double x2 = TheNet.Node[Edge[1]][0];
			double y2 = TheNet.Node[Edge[1]][1];
			double hm = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
			Matrix A = {
				{2,1},
				{1,2}
			};

			double mult = Betta((int)Edge[3]) * hm / 6.;

			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					A[i][j] = A[i][j] * mult;
				}
			}

			vector<double> b = {
					mult * (2 * UB(TheNet.Node[Edge[0]],Edge[2]) + UB(TheNet.Node[Edge[1]],Edge[2])),
					mult * (UB(TheNet.Node[Edge[0]],Edge[2]) + 2 * UB(TheNet.Node[Edge[1]],Edge[2]))
			};
			ToGLobalProf(A, b, Edge);
			ToGlobalPlot(A, b, Edge);
		}
	}

	void AddSecondCondi()
	{
		int length = TheNet.SecondCondi.size();
		for (int i = 0; i < length; i++)
		{
			vector<int> Edge = TheNet.SecondCondi[i];
			double x1 = TheNet.Node[Edge[0]][0];
			double y1 = TheNet.Node[Edge[0]][1];
			double x2 = TheNet.Node[Edge[1]][0];
			double y2 = TheNet.Node[Edge[1]][1];
			double hm = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

			double mult = hm / 6.;
			double toDelete = Tetta(TheNet.Node[Edge[0]], Edge[2]);
			vector<double> b = {
					mult * (2 * Tetta(TheNet.Node[Edge[0]],Edge[2]) + Tetta(TheNet.Node[Edge[1]],Edge[2])),
					mult * (Tetta(TheNet.Node[Edge[0]],Edge[2]) + 2 * Tetta(TheNet.Node[Edge[1]],Edge[2]))
			};
			this->b[Edge[0]] += b[0];
			this->b[Edge[1]] += b[1];
		}
	}

	void AddFirst()
	{
		double max = 0;
		int length = AProf.AL.size();
		for (int i = 0; i < length; i++)
		{
			if (max < abs(AProf.AL[i]))
			{
				max = abs(AProf.AL[i]);
			}
		}
		max *= 1e+30;

		length = TheNet.firstCondi.size();
		for (int i = 0; i < length; i++)
		{
			int n = TheNet.firstCondi[i][0];
			AProf.DI[n] = max;
			A[n][n] = max; //ToRemove
			b[n] = max * Ug(TheNet.Node[n], TheNet.firstCondi[i][1]);
			q[n] = Ug(TheNet.Node[n], TheNet.firstCondi[i][1]);
		}
	}
	//Решение Ситемы

	void LUFactorization(MatrixProf& A, MatrixProf& LU)
	{
		int length = A.IA.size();
		for (int i = 0; i < length; i++)
		{
			A.IA[i]--;
		}
		LU.size = A.size;
		LU.IA.resize(LU.size + 1);

		for (int i = 0; i < A.size + 1; i++)
			LU.IA[i] = A.IA[i];

		LU.AL.resize(LU.IA[LU.size]);
		LU.AU.resize(LU.IA[LU.size]);
		LU.JA.resize(LU.IA[LU.size]);
		LU.DI.resize(LU.size);

		for (int i = 0; i < A.IA[A.size]; i++)
			LU.JA[i] = A.JA[i];

		for (int i = 0; i < A.size; i++)
		{
			double sumD = 0;
			int i0 = A.IA[i], i1 = A.IA[i + 1];

			for (int k = i0; k < i1; k++)
			{
				double sumL = 0, sumU = 0;
				int j = A.JA[k];

				// Calculate L[i][j], U[j][i]
				int j0 = A.IA[j], j1 = A.IA[j + 1];

				int kl = i0, ku = j0;

				for (; kl < i1 && ku < j1; )
				{
					int j_kl = A.JA[kl];
					int j_ku = A.JA[ku];

					if (j_kl == j_ku)
					{
						sumL += LU.AL[kl] * LU.AU[ku];
						sumU += LU.AU[kl] * LU.AL[ku];
						kl++;
						ku++;
					}
					if (j_kl > j_ku)
						ku++;
					if (j_kl < j_ku)
						kl++;
				}

				LU.AL[k] = A.AL[k] - sumL;
				LU.AU[k] = A.AU[k] - sumU;
				LU.AU[k] /= A.DI[j];

				// Calculate sum for DI[i]
				sumD += LU.AL[k] * LU.AU[k];
			}

			// Calculate DI[i]
			LU.DI[i] = A.DI[i] - sumD;
		}
	}

	Matrix ProfToPlot(MatrixProf& A)
	{
		Matrix Res(A.size);
		int n = A.size;
		for (int i = 0; i < n; i++)
		{
			Res[i].resize(n);
			Res[i][i] = A.DI[i];
		}

		for (int i = 0; i < n; i++)
		{
			for (int jadr = A.IA[i]; jadr < A.IA[i + 1]; jadr++)
			{
				int j = A.JA[jadr];
				Res[i][j] = A.AL[jadr];
				Res[j][i] = A.AU[jadr];
			}
		}
		return Res;
	}

	void Calculate()
	{
		LUFactorization(AProf, LU);
		AuxVectors TmpSolution;
		TmpSolution.Ax = vector<double>(AProf.size);
		TmpSolution.LU = vector<double>(AProf.size);
		TmpSolution.p = vector<double>(AProf.size);
		TmpSolution.r = vector<double>(AProf.size);
		TmpSolution.z = vector<double>(AProf.size);
		TmpSolution.temp = vector<double>(AProf.size);
		LOS_LU(AProf, q, b, LU, TmpSolution, 10000, 1e-13);
	}

	void LOS_LU(MatrixProf& A, vector<double>& x, vector<double>& f, MatrixProf& LU, AuxVectors& aux, int maxiter, double eps)
	{
		int size = A.size;

		// Calculate r0
		Multiply(A, x, aux.Ax);
		for (int i = 0; i < size; i++)
			aux.r[i] = f[i] - aux.Ax[i];
		Forward(LU, aux.r, aux.r);

		//Calculate z0
		Backward(LU, aux.z, aux.r);

		// Calculate p0
		Multiply(A, aux.z, aux.p);
		Forward(LU, aux.p, aux.p);

		double diff = MultVecs(size, aux.r, aux.r);

		int k = 0;
		for (; k < maxiter && diff >= eps; k++)
		{
			// Calculate alpha
			double dotP = MultVecs(size, aux.p, aux.p);
			double a = MultVecs(size, aux.p, aux.r) / dotP;

			// Calculate xk, rk
			for (int i = 0; i < size; i++)
			{
				x[i] += a * aux.z[i];
				aux.r[i] -= a * aux.p[i];
			}

			// Calculate beta
			Backward(LU, aux.Ax, aux.r);
			Multiply(A, aux.Ax, aux.temp);
			Forward(LU, aux.Ax, aux.temp);
			double b = -MultVecs(size, aux.p, aux.Ax) / dotP;

			// Calculate zk, pk
			Backward(LU, aux.temp, aux.r);
			for (int i = 0; i < size; i++)
			{
				aux.z[i] = aux.temp[i] + b * aux.z[i];
				aux.p[i] = aux.Ax[i] + b * aux.p[i];
			}

			// Calculate difference
			diff = MultVecs(size, aux.r, aux.r);
		}
		maxiter = k;
	}

private:

	double Ug(vector<double>& node, int k)
	{
		double x = node[0];
		double y = node[1];
		return x + 1;
	}

	double UB(vector<double>& node, int k)
	{
		double x = node[0];
		double y = node[1];
		return -1.5 * y * y + y;
	}

	double Tetta(vector<double>& node, int k)
	{
		double x = node[0];
		double y = node[1];
		return 2 * x * x + 3 * y * y;
	}


	double F(double x, double y, int field)
	{
		return -x - 3 * y;
	}

	double Lambda(vector<double>& node, int field)
	{
		double x = node[0];
		double y = node[1];
		return 2 * x * x + 3 * y * y;
	}

	double Betta(int field)
	{
		return 2;
	}

	double Gamma(int field)
	{
		return 3;
	}
	//параметры

	void ToGlobalPlot(Matrix& L, vector<double>& b, vector<int>& el)
	{
		int length = L.size();
		for (int i = 0; i < length; i++)
		{
			for (int j = 0; j < length; j++)
			{
				A[el[i]][el[j]] += L[i][j];
			}
		}

		//for (int i = 0; i < length; i++)
		//{
		//	this->b[el[i]] += b[i];
		//}
	}

	void ToGLobalProf(Matrix& A, vector<double>& b, vector<int>& el)
	{
		int length = A.size();
		for (int i = 0; i < length; i++)
		{
			AProf.DI[el[i]] = AProf.DI[el[i]] + A[i][i];
		}

		for (int i = 0; i < length; i++)
		{
			int ibeg = AProf.IA[el[i]] - 1;
			for (int j = 0; j < i; j++)
			{
				int iend = AProf.IA[el[i] + 1] - 1;
				while (AProf.JA[ibeg] != el[j])
				{
					int ind = (ibeg + iend) / 2;
					if (AProf.JA[ind] <= el[j])
					{
						ibeg = ind;
					}
					else
					{
						iend = ind;
					}
				}
				AProf.AL[ibeg] += A[i][j];
				AProf.AU[ibeg] += A[j][i];
				ibeg++;
			}
		}

		for (int i = 0; i < length; i++)
		{
			this->b[el[i]] += b[i];
		}

	}

	double MultVecs(int size, vector<double>& vec1, vector<double>& vec2)
	{
		double sum = 0;
		for (int i = 0; i < size; i++)
			sum += vec1[i] * vec2[i];
		return sum;
	}

	void Multiply(MatrixProf& A, vector<double>& vec, vector<double>& res)
	{
		int size = A.size;

		for (int i = 0; i < size; i++)
		{
			res[i] = vec[i] * A.DI[i];
			for (int k = A.IA[i]; k < A.IA[i + 1]; k++)
			{
				int j = A.JA[k];
				res[i] += A.AL[k] * vec[j];
				res[j] += A.AU[k] * vec[i];
			}
		}
	}

	void Forward(MatrixProf& A, vector<double>& x, vector<double>& b)
	{
		int size = A.size;

		for (int i = 0; i < size; i++)
		{
			double sum = 0;
			int i0 = A.IA[i], i1 = A.IA[i + 1];
			for (int k = i0; k < i1; k++)
			{
				int j = A.JA[k];
				sum += A.AL[k] * x[j];
			}
			x[i] = (b[i] - sum) / A.DI[i];
		}
	}

	void Backward(MatrixProf& A, vector<double>& x, vector<double>& b)
	{
		int size = A.size;

		for (int i = 0; i < size; i++)
			x[i] = b[i];

		for (int i = size - 1; i >= 0; i--)
		{
			int i0 = A.IA[i], i1 = A.IA[i + 1];
			for (int k = i0; k < i1; k++)
			{
				int j = A.JA[k];
				x[j] -= A.AU[k] * x[i];
			}
		}
	}
};








int main()
{
	//Nett.BuildNet(0, 2, 0, 2, 2, 2);
	//Nett.Node = vector<vector<double>>{
	//{2,0},
	//{2,1},
	//{3,1},
	//{2,4},
	//{7,4}
	//};
	//Nett.Elements = vector<vector<int>>{
	//	{0,1,2},
	//	{2,3,4},
	//	{1,2,3}
	//};
	//Nett.fields = vector<int>(3);
	//Nett.fields[1] = 1;
	//Nett.fields[2] = 1;
	//Nett.ThirdCondi.resize(1);
	//Nett.ThirdCondi[0] = { 2,4,0,0};
	//Nett.SecondCondi = {
	//	///*{0,1,1},
	//	//*/{1,3,1},
	//	{3,4,0}
	//};
	//Nett.firstCondi = {
	//	{0,0},
	//	{2,0}
	//};

	fstream nodes;
	fstream elements;
	fstream fields;
	fstream condi1;
	fstream condi2;
	fstream condi3;
	nodes.open("nodes.txt");
	elements.open("elements.txt");
	fields.open("fields.txt");
	condi1.open("condi1.txt");
	condi2.open("condi2.txt");
	condi3.open("condi3.txt");
	
	//Net Nett;
	//Nett.BuildNet(1, 2, 1, 2, 2, 2);
	//Nett.SaveNet(nodes, elements, fields);
	//nodes.close();
	//elements.close();
	//fields.close();
	Net Nett(nodes,elements,fields, condi1,condi2,condi3);

	//Nett.Node = { {1,1}, {2,1}, {1.5,1.5}, {1,2},{2,2} };
	//Nett.Elements = { {0,1,2},{0,2,3},{1,2,4},{2,3,4} };
	//Nett.firstCondi = { {0,0},{1,0} };
	//Nett.SecondCondi = {
	//	{1,4,0},
	//	{3,4,1},
	//	//{0,3,2}
	//};
	//Nett.ThirdCondi = {
	//	{0,3,0,0}
	//};
	//Nett.fields.resize(Nett.Elements.size());
	//for (size_t i = 0; i < Nett.Elements.size(); i++)
	//{
	//	Nett.fields[i] = 0;
	//}


	cout << scientific << setprecision(15);
	Eq Equation(Nett);
	Equation.BuildGlobalProf();
	Equation.AddThirdCondi();
	Equation.AddSecondCondi();
	Equation.AddFirst();
	
	/*Equation.LUFactorization(Equation.AProf, Equation.LU);
	AuxVectors TmpSolution;
	TmpSolution.Ax = vector<double>(Equation.AProf.size);
	TmpSolution.LU = vector<double>(Equation.AProf.size);
	TmpSolution.p = vector<double>(Equation.AProf.size);
	TmpSolution.r = vector<double>(Equation.AProf.size);
	TmpSolution.z = vector<double>(Equation.AProf.size);
	TmpSolution.temp = vector<double>(Equation.AProf.size);
	Equation.LOS_LU(Equation.AProf, Equation.q, Equation.b, Equation.LU, TmpSolution, 10000, 1e-13);*/
	//Equation.PrintPlot(Equation.A);

	Equation.Calculate();
	for (size_t i = 0; i < Equation.AProf.size; i++)
	{
		double U = Equation.U(Equation.TheNet.Node[i][0], Equation.TheNet.Node[i][1]);
		cout << Equation.q[i] << " " << U << " " << abs(U - Equation.q[i]) << endl;
	}

	std::cout << "Hello World!\n";
}