const int Node_Num = 9;
const int Element_Num = 8;
const int Global_m = 1;

const double PI = 3.1415926;
const double Gauss_T = 0.5773502692;
const double Gauss_T2 = 0.7745966692;
const double Grav = 9.8;
const double Omega = 2 * PI;

enum Boud_Type { NONE, FREE, RIGID, CORNER };

struct Node_Type
{
	double x, y;
	Boud_Type boudary;
};

struct Element_Type
{
	unsigned int p[3];
};

Node_Type node[Node_Num] =
{
{1, 0, FREE},  {1, 1, FREE},  {1, 2, FREE},
{2, 0, RIGID}, {2, 1, NONE},  {2, 2, RIGID},
{3, 0, RIGID}, {3, 1, RIGID}, {3, 2, RIGID}
};

Element_Type element[Element_Num] =
{
{1, 4, 2}, {4, 5, 2}, {4, 7, 5}, {7, 8, 5},
{2, 5, 3}, {5, 6, 3}, {5, 8, 6}, {8, 9, 6}
};

unsigned int Boudary_Num;

double Ae[3][3], Be[3][3], Ce[3][3], De[3][3];

double A[Node_Num][Node_Num],
       B[Node_Num][Node_Num],
       C[Node_Num][Node_Num],
       D[Node_Num][Node_Num];

double G[3 * Node_Num][3 * Node_Num];

void boud_elem(void)
{
	Boudary_Num = 0;
	for(int i = 0; i < Node_Num; i++) Boudary_Num += (node[i].boudary == FREE) ? 1 : 0;
}

void gen_elem(unsigned int el)
{
	int i, j, k;
	double gauss;
	double dx, dy;
	double swap, temp;
	double er, ez, e_k, e_b;

	double ea[3], eb[3], ec[3];
	double x[3], y[3];
	double sx[3], sy[3];
	double nt_x[3], nt_y[3];
	double xc, yc;
	double nr[3];
	double s;
	double edge_len[3];
	Boud_Type bt[3];
	Node_Type *np[3];
//	Node_Type *nsp[3];

	for(i = 0; i < 3; i++)	//Get x, y of all three point.
	{
		np[i] = &node[element[el].p[i]];
		x[i] = sx[i] = np[i]->x;
		y[i] = sy[i] = np[i]->y;
	}

	for(i = 0; i < 3; i++)	//Sorting the point on x.
	{
		for(j = 0; j < 2; j++)
		{
			if (sx[j] > sx[j+1])
			{
				swap = sx[j];
				sx[j] = sx[j+1];
				sx[j+1] = swap;
			}
		}
	}

	for(i = 0; i < 3; i++)
	{
		nt_x[i] = x[(i+1)%3];
		nt_y[i] = y[(i+1)%3];
		dx = nt_x[i] - x[i];
		dy = nt_y[i] - y[i];
		edge_len[i] = sqrt(dx * dx + dy * dy);

		if ((np[i]->boudary == NONE)||(np[(i+1)%3]->boudary == NONE)) bt[i] = NONE;
		else
		{
			if ((np[i]->boudary == FREE)&&(np[(i+1)%3]->boudary == FREE)) bt[i] = FREE;
			else bt[i] = RIGID;
		}
		nr[i] = dy / edge_len[i];
	}

	xc = (x[0] + x[1] + x[2]) / 3.0;
	yc = (y[0] + y[1] + y[2]) / 3.0;

	ea[0] = x[1] * y[2] - x[2] * y[1];
	ea[1] = x[2] * y[0] - x[0] * y[2];
	ea[2] = x[0] * y[1] - x[1] * y[0];

	eb[0] = y[1] - y[2];
	eb[1] = y[2] - y[0];
	eb[2] = y[0] - y[1];

	ec[0] = x[2] - x[1];
	ec[1] = x[0] - x[2];
	ec[2] = x[1] - x[0];

	s = 0.5 * (ea[0] + ea[1] + ea[2]);


	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			double ez1, ez2;
			double eA, eB, eC;

			swap = temp = 0;

			if (sx[1] != sx[0])
			{
				swap = temp = 0;

				e_k = 2 / (sx[1] - sx[0]);
				e_b = (sx[0] + sx[1]) / (sx[0] - sx[1]);

				gauss = Gauss_T;
				er = (gauss - e_b) / e_k;
				ez1 = (sy[2] - sy[0]) * (er - sx[0]) / (sx[2] - sx[0]) + sy[0];
				ez2 = (sy[1] - sy[0]) * (er - sx[0]) / (sx[1] - sx[0]) + sy[0];
				eA = (ea[i] + eb[i] * er) * (ea[j] + eb[j] * er);
				eB = ((ea[i] + eb[i] * er) * ec[j] + (ea[j] + eb[j] * er) * ec[i]) / 2;
				eC = (ec[i] * ec[j]) / 3;
				temp = (eA * fabs(ez2 - ez1) + eB * fabs(ez2 * ez2 - ez1 * ez1) + eC * fabs(ez2 * ez2 * ez2 - ez1 * ez1 * ez1)) / er / er;

				gauss = -1 * Gauss_T;
				er = (gauss - e_b) / e_k;
				ez1 = (sy[2] - sy[0]) * (er - sx[0]) / (sx[2] - sx[0]) + sy[0];
				ez2 = (sy[1] - sy[0]) * (er - sx[0]) / (sx[1] - sx[0]) + sy[0];
				eA = (ea[i] + eb[i] * er) * (ea[j] + eb[j] * er);
				eB = ((ea[i] + eb[i] * er) * ec[j] + (ea[j] + eb[j] * er) * ec[i]) / 2;
				eC = (ec[i] * ec[j]) / 3;
				temp += (eA * fabs(ez2 - ez1) + eB * fabs(ez2 * ez2 - ez1 * ez1) + eC * fabs(ez2 * ez2 * ez2 - ez1 * ez1 * ez1)) / er / er;

				swap = temp * Global_m * Global_m * (sx[1] - sx[0]) / 8 / s / s;

				Ae[i][j] += swap;
			}

			if (sx[1] != sx[2])
			{
				swap = temp = 0;

				e_k = 2 / (sx[1] - sx[2]);
				e_b = (sx[2] + sx[1]) / (sx[2] - sx[1]);

				gauss = Gauss_T;
				er = (gauss - e_b) / e_k;
				ez1 = (sy[2] - sy[0]) * (er - sx[0]) / (sx[2] - sx[0]) + sy[0];
				ez2 = (sy[1] - sy[2]) * (er - sx[2]) / (sx[1] - sx[2]) + sy[2];
				eA = (ea[i] + eb[i] * er) * (ea[j] + eb[j] * er);
				eB = ((ea[i] + eb[i] * er) * ec[j] + (ea[j] + eb[j] * er) * ec[i]) / 2;
				eC = (ec[i] * ec[j]) / 3;
				temp = (eA * fabs(ez2 - ez1) + eB * fabs(ez2 * ez2 - ez1 * ez1) + eC * fabs(ez2 * ez2 * ez2 - ez1 * ez1 * ez1)) / er / er;

				gauss = -1 * Gauss_T;
				er = (gauss - e_b) / e_k;
				ez1 = (sy[2] - sy[0]) * (er - sx[0]) / (sx[2] - sx[0]) + sy[0];
				ez2 = (sy[1] - sy[2]) * (er - sx[2]) / (sx[1] - sx[2]) + sy[2];
				eA = (ea[i] + eb[i] * er) * (ea[j] + eb[j] * er);
				eB = ((ea[i] + eb[i] * er) * ec[j] + (ea[j] + eb[j] * er) * ec[i]) / 2;
				eC = (ec[i] * ec[j]) / 3;
				temp += (eA * fabs(ez2 - ez1) + eB * fabs(ez2 * ez2 - ez1 * ez1) + eC * fabs(ez2 * ez2 * ez2 - ez1 * ez1 * ez1)) / er / er;

				swap = temp * Global_m * Global_m * (sx[1] - sx[2]) / 8 / s / s;

				Ae[i][j] += swap;
			}

			Ce[i][j] = ec[i] * ec[j] * xc / 4.0 / s;
			Ae[i][j] = eb[i] * eb[j] * xc / 4.0 / s + Ce[i][j];
		}

		if ((bt[i] == FREE)||(bt[i] == RIGID))
		{
			if (nt_x[i] != x[i])
			{
				e_k = 2 / (nt_x[i] - x[i]);
				e_b = (nt_x[i] + x[i]) / (x[i] - nt_x[i]);

				gauss = Gauss_T;
				er = (gauss - e_b) / e_k;
				ez = (nt_y[i] - y[i]) * (er - x[i]) / (nt_x[i] - x[i]) + y[i];
				Be[i][(i+1)%3] = (ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez);

				gauss = -1 * Gauss_T;
				er = (gauss - e_b) / e_k;
				ez = (nt_y[i] - y[i]) * (er - x[i]) / (nt_x[i] - x[i]) + y[i];
				Be[i][(i+1)%3] += (ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez);

				Be[i][(i+1)%3] *= nr[i] * edge_len[i] / 8 / s / s;
			}
			else
			{
				e_k = 2 / (nt_y[i] - y[i]);
				e_b = (nt_y[i] + y[i]) / (y[i] - nt_y[i]);

				gauss = Gauss_T;
				ez = (gauss - e_b) / e_k;
				er = (nt_x[i] - x[i]) * (ez - y[i]) / (nt_y[i] - y[i]) + x[i];
				Be[i][(i+1)%3] = (ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez);

				gauss = -1 * Gauss_T;
				ez = (gauss - e_b) / e_k;
				er = (nt_x[i] - x[i]) * (ez - y[i]) / (nt_y[i] - y[i]) + x[i];
				Be[i][(i+1)%3] += (ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez);

				Be[i][(i+1)%3] *= nr[i] * edge_len[i] / 8 / s / s;

			}
		}

		if (bt[i] == FREE)
		{
			if (nt_x[i] != x[i])
			{
				e_k = 2 / (nt_x[i] - x[i]);
				e_b = (nt_x[i] + x[i]) / (x[i] - nt_x[i]);

				gauss = Gauss_T;
				er = (gauss - e_b) / e_k;
				ez = (nt_y[i] - y[i]) * (er - x[i]) / (nt_x[i] - x[i]) + y[i];
				De[i][(i+1)%3] = ((ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez)) / sqrt(1 + (Grav / er / Omega / Omega) * (Grav / er / Omega / Omega));

				gauss = -1 * Gauss_T;
				er = (gauss - e_b) / e_k;
				ez = (nt_y[i] - y[i]) * (er - x[i]) / (nt_x[i] - x[i]) + y[i];
				De[i][(i+1)%3] += ((ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez)) / sqrt(1 + (Grav / er / Omega / Omega) * (Grav / er / Omega / Omega));

				De[i][(i+1)%3] *= edge_len[i] / 8 / s / s;
			}
			else
			{
				e_k = 2 / (nt_y[i] - y[i]);
				e_b = (nt_y[i] + y[i]) / (y[i] - nt_y[i]);

				gauss = Gauss_T;
				ez = (gauss - e_b) / e_k;
				er = (nt_x[i] - x[i]) * (ez - y[i]) / (nt_y[i] - y[i]) + x[i];
				De[i][(i+1)%3] = ((ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez)) / sqrt(1 + (Grav / er / Omega / Omega) * (Grav / er / Omega / Omega));

				gauss = -1 * Gauss_T;
				ez = (gauss - e_b) / e_k;
				er = (nt_x[i] - x[i]) * (ez - y[i]) / (nt_y[i] - y[i]) + x[i];
				De[i][(i+1)%3] += ((ea[i] + eb[i] * er + ec[i] * ez) * (ea[(i+1)%3] + eb[(i+1)%3] * er + ec[(i+1)%3] * ez)) / sqrt(1 + (Grav / er / Omega / Omega) * (Grav / er / Omega / Omega));

				De[i][(i+1)%3] *= edge_len[i] / 8 / s / s;

			}

		}
	}
}

void gen_L(void)
{

}

void gen_total(void)
{

}