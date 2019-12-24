#include <iostream>
#include <vector>
#include <functional>
#include <iomanip>
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
		//разбиение на подобласти
		//DevideBy2Fields();
	}
private:
};


class Eq
{
public:

	Net FuckingNet;
	double lambda = 1, gamma = 0;
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
		return 2 * x * x + 3 * y * y;
	}

	double Ug(vector<double>& node, int k)
	{
		double x = node[0];
		double y = node[1];
		return 2 * x * x + 3;
	}

	double UB(vector<double>& node, int k)
	{
		double x = node[0];
		double y = node[1];
		return -4 * x + 2 + 3 * y * y;
	}

	double Tetta(vector<double>& node, int k)
	{
		double x = node[0];
		double y = node[1];
		switch (k)
		{
		case 0:
			return 8*x;
		case 1:
			return 12*y;
		default:
			break;
		}
	}


	double F(double x, double y, int field)
	{
		return 6 * x*x + 9 * y*y -20;
	}

	double Lambda(vector<double>& node, int field)
	{
		return 2;
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

	Eq()
	{
		FuckingNet.BuildNet(0, 2, 0, 2, 2, 2);
		A = Matrix(FuckingNet.Node.size());
		b = vector<double>(FuckingNet.Node.size());
	}

	Eq(Net net)
	{
		FuckingNet = net;
		A = Matrix(FuckingNet.Node.size());
		for (int i = 0; i < A.size(); i++)
		{
			A[i] = vector<double>(A.size());
		}
		b = vector<double>(FuckingNet.Node.size());
		q = vector<double>(FuckingNet.Node.size());
	}
	//построение матрицы G
	vector<vector<double>> BuildG(vector<vector<double>>& D_1, double DetD, vector<int>& el, int field) {
		vector<vector<double>> G(3);
		double multix = abs(DetD) / 2;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{

				G[i].push_back(Lambda(FuckingNet.Node[0], field) * multix * (D_1[i][1] * D_1[j][1] + D_1[i][2] * D_1[j][2])); // Lambda = const;
			}
		}
		return G;
	}
	//разложение коэф дифузии по квадратичным базисным функциям
	vector<vector<double>> BuildGDecomposeSquardL(vector<vector<double>>& D_1, double DetD, vector<int>& el, int field) {
		vector<vector<double>> G(3);
		double multix = abs(DetD) / 6;
		double x1 = FuckingNet.Node[el[0]][0];
		double y1 = FuckingNet.Node[el[0]][1];
		double x2 = FuckingNet.Node[el[1]][0];
		double y2 = FuckingNet.Node[el[1]][1];
		double x3 = FuckingNet.Node[el[2]][0];
		double y3 = FuckingNet.Node[el[2]][1];
		vector<double> point4 = { findMax(x1, x2) - abs(x1 - x2) / 1, findMax(y1,y2) - abs(y1 - y2) };
		vector<double> point5 = { findMax(x3, x2) - abs(x3 - x2) / 1, findMax(y3,y2) - abs(y3 - y2) };
		vector<double> point6 = { findMax(x1, x3) - abs(x1 - x3) / 1, findMax(y1,y3) - abs(y1 - y3) };

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double sumL = Lambda(point4, field) +
					Lambda(point5, field) +
					Lambda(point6, field);
				G[i].push_back(sumL * multix * (D_1[i][1] * D_1[j][1] + D_1[i][2] * D_1[j][2])); // Lambda = const;
			}
		}
		return G;
	}
	//разложение коэф дифузии по линейным базисным функциям
	vector<vector<double>> BuildGDecomposeLinalL(vector<vector<double>>& D_1, double DetD, vector<int>& el, int field) {
		vector<vector<double>> G(3);
		double multix = abs(DetD) / 6;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				double sumL = Lambda(FuckingNet.Node[el[0]], field) +
					Lambda(FuckingNet.Node[el[1]], field) +
					Lambda(FuckingNet.Node[el[2]], field);
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

	void BuildGlobalProf()
	{
		BuildProfile();

		for (int i = 0; i < FuckingNet.Elements.size(); i++)
		{
			vector<int> element = FuckingNet.Elements[i];
			int field = FuckingNet.fields[i];
			BuildLocal(element, field);
			PrintPlot(A);
			cout << "\n\n\n";
		}
	}

	Matrix BuildLocal(vector<int>& el, int field)
	{
		double x1 = FuckingNet.Node[el[0]][0];
		double x2 = FuckingNet.Node[el[1]][0];
		double x3 = FuckingNet.Node[el[2]][0];
		double y1 = FuckingNet.Node[el[0]][1];
		double y2 = FuckingNet.Node[el[1]][1];
		double y3 = FuckingNet.Node[el[2]][1];

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
		vector<vector<double>> G = BuildGDecomposeLinalL(D_1, DetD, el, field);
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
		for (int i = 0; i < FuckingNet.Elements.size(); i++)
		{
			Matrix test = BuildLocal(FuckingNet.Elements[i], FuckingNet.fields[i]);
		}
	}

	void BuildProfile()
	{
		vector<vector<int>> profile(FuckingNet.Node.size());

		for (int i = 0; i < FuckingNet.Elements.size(); i++)
		{
			for (int j = 1; j < 3; j++)
			{
				for (int k = 0; k < j; k++)
				{
					int current = FuckingNet.Elements[i][j];
					int node = FuckingNet.Elements[i][k];
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
		for (int i = 1; i < FuckingNet.Node.size(); i++)
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
		AProf.DI = vector<double>(FuckingNet.Node.size());
		AProf.size = AProf.DI.size();
	}

	void AddThirdCondi()
	{
		int length = FuckingNet.ThirdCondi.size();
		for (int i = 0; i < length; i++)
		{
			vector<int> Edge = FuckingNet.ThirdCondi[i];
			double x1 = FuckingNet.Node[Edge[0]][0];
			double y1 = FuckingNet.Node[Edge[0]][1];
			double x2 = FuckingNet.Node[Edge[1]][0];
			double y2 = FuckingNet.Node[Edge[1]][1];
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
					mult * (2 * UB(FuckingNet.Node[Edge[0]],Edge[2]) + UB(FuckingNet.Node[Edge[1]],Edge[2])),
					mult * (UB(FuckingNet.Node[Edge[0]],Edge[2]) + 2 * UB(FuckingNet.Node[Edge[1]],Edge[2]))
			};
			ToGLobalProf(A, b, Edge);
			ToGlobalPlot(A, b, Edge);
		}
	}

	void AddSecondCondi()
	{
		int length = FuckingNet.SecondCondi.size();
		for (int i = 0; i < length; i++)
		{
			vector<int> Edge = FuckingNet.SecondCondi[i];
			double x1 = FuckingNet.Node[Edge[0]][0];
			double y1 = FuckingNet.Node[Edge[0]][1];
			double x2 = FuckingNet.Node[Edge[1]][0];
			double y2 = FuckingNet.Node[Edge[1]][1];
			double hm = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));

			double mult = hm / 6.;
			double toDelete = Tetta(FuckingNet.Node[Edge[0]], Edge[2]);
			vector<double> b = {
					mult * (2 * Tetta(FuckingNet.Node[Edge[0]],Edge[2]) + Tetta(FuckingNet.Node[Edge[1]],Edge[2])),
					mult * (Tetta(FuckingNet.Node[Edge[0]],Edge[2]) + 2 * Tetta(FuckingNet.Node[Edge[1]],Edge[2]))
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

		length = FuckingNet.firstCondi.size();
		for (int i = 0; i < length; i++)
		{
			int n = FuckingNet.firstCondi[i][0];
			AProf.DI[n] = max;
			A[n][n] = max; //ToRemove
			b[n] = max * Ug(FuckingNet.Node[n], FuckingNet.firstCondi[i][1]);
			q[n] = Ug(FuckingNet.Node[n], FuckingNet.firstCondi[i][1]);
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
};








int main()
{
	Net testNet;
	//testNet.BuildNet(0, 2, 0, 2, 2, 2);
	//testNet.Node = vector<vector<double>>{
	//{2,0},
	//{2,1},
	//{3,1},
	//{2,4},
	//{7,4}
	//};
	//testNet.Elements = vector<vector<int>>{
	//	{0,1,2},
	//	{2,3,4},
	//	{1,2,3}
	//};
	//testNet.fields = vector<int>(3);
	//testNet.fields[1] = 1;
	//testNet.fields[2] = 1;
	//testNet.ThirdCondi.resize(1);
	//testNet.ThirdCondi[0] = { 2,4,0,0};
	//testNet.SecondCondi = {
	//	///*{0,1,1},
	//	//*/{1,3,1},
	//	{3,4,0}
	//};
	//testNet.firstCondi = {
	//	{0,0},
	//	{2,0}
	//};
	//testNet.BuildNet(0,2,0,2,2,2);


	testNet.Node = { {1,1}, {2,1}, {1.5,1.5}, {1,2},{2,2} };
	testNet.Elements = { {0,1,2},{0,2,3},{1,2,4},{2,3,4} };
	testNet.firstCondi = { {0,0},{1,0} };
	testNet.SecondCondi = {
		{1,4,0},
		{3,4,1},
		//{0,3,2}
	};
	testNet.ThirdCondi = {
		{0,3,0,0}
	};
	testNet.fields.resize(testNet.Elements.size());
	for (size_t i = 0; i < testNet.Elements.size(); i++)
	{
		testNet.fields[i] = 0;
	}


	cout << scientific << setprecision(15);
	Eq testSys(testNet);
	testSys.BuildGlobalProf();
	testSys.AddThirdCondi();
	testSys.AddSecondCondi();
	testSys.AddFirst();


	testSys.PrintPlot(testSys.A);
	testSys.LUFactorization(testSys.AProf, testSys.LU);
	AuxVectors useless;
	useless.Ax = vector<double>(testSys.AProf.size);
	useless.LU = vector<double>(testSys.AProf.size);
	useless.p = vector<double>(testSys.AProf.size);
	useless.r = vector<double>(testSys.AProf.size);
	useless.z = vector<double>(testSys.AProf.size);
	useless.temp = vector<double>(testSys.AProf.size);
	testSys.LOS_LU(testSys.AProf, testSys.q, testSys.b, testSys.LU, useless, 10000, 1e-13);
	for (size_t i = 0; i < testSys.AProf.size; i++)
	{
		double U = testSys.U(testSys.FuckingNet.Node[i][0], testSys.FuckingNet.Node[i][1]);
		cout << testSys.q[i] << " " << U << " " << abs(U - testSys.q[i]) << endl;
	}
	vector<double> rezult(testSys.AProf.size);
	testSys.Multiply(testSys.AProf, testSys.q, rezult);
	std::cout << "Hello World!\n";
}