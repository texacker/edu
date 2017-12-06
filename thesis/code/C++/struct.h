const int Node_Num = 30;		//const int Edge_Num = 16;
const int Element_Num = 32;
const int Free_Num = 10;

const double PI = 3.1415926;
const double Gauss_T = 0.5773502692;
const double Gauss_T2 = 0.7745966692;
const double Grav = 9.8;

const int Global_m = 1;
const double Eps = 0.0000001;
const double Omega = 2 * PI;
const double Fea_Len = 1;
const double Gause = Grav / Omega / Omega / Fea_Len;

double Ae[3][3], Be[3][3], Ce[3][3], De[3][3];

double A[Node_Num][Node_Num],
       B[Node_Num][Node_Num],
       C[Node_Num][Node_Num],
       D[Node_Num][Node_Num];

double vr[2 * (Node_Num + Free_Num)];
double vi[2 * (Node_Num + Free_Num)];
double far G[2 * (Node_Num + Free_Num)][2 * (Node_Num + Free_Num)];

//enum Boud_Type { NONE, FREE, RIGID, CORNER };
enum Edge_Type { INSIDE, FREE, RIGID };

struct Node_Type
{
	double x, y;
};

struct Element_Type
{
	unsigned int n[3];
	Edge_Type e[3];
};

struct Element_Data
{
	double x[3], y[3];
	double sx[3], sy[3];
	double nx[3], ny[3];
	double dx[3], dy[3];
	double s, xc, yc;
	double len[3], nr[3];
	double a[3], b[3], c[3];
};

Node_Type node[Node_Num] =
{
{ 2, 0}, { 2, 0.6}, { 2, 1.2}, { 2, 1.8}, { 2, 2.4},
{-2, 0}, {-2, 0.6}, {-2, 1.2}, {-2, 1.8}, {-2, 2.4},
{ 3, 0}, { 3, 0.6}, { 3, 1.2}, { 3, 1.8}, { 3, 2.4},
{ 4, 0}, { 4, 0.6}, { 4, 1.2}, { 4, 1.8}, { 4, 2.4},
{-3, 0}, {-3, 0.6}, {-3, 1.2}, {-3, 1.8}, {-3, 2.4},
{-4, 0}, {-4, 0.6}, {-4, 1.2}, {-4, 1.8}, {-4, 2.4}
};

Element_Type element[Element_Num] =
{
{{ 0, 10, 11}, {RIGID, INSIDE, INSIDE}},
{{ 0, 11,  1}, {INSIDE, INSIDE, FREE}},
{{ 1, 11, 12}, {INSIDE, INSIDE, INSIDE}},
{{ 1, 12,  2}, {INSIDE, INSIDE, FREE}},
{{ 2, 12, 13}, {INSIDE, INSIDE, INSIDE}},
{{ 2, 13,  3}, {INSIDE, INSIDE, FREE}},
{{ 3, 13, 14}, {INSIDE, INSIDE, INSIDE}},
{{ 3, 14,  4}, {INSIDE, RIGID, FREE}},

{{10, 15, 16}, {RIGID, RIGID, INSIDE}},
{{10, 16, 11}, {INSIDE, INSIDE, INSIDE}},
{{11, 16, 17}, {INSIDE, RIGID, INSIDE}},
{{11, 17, 12}, {INSIDE, INSIDE, INSIDE}},
{{12, 17, 18}, {INSIDE, RIGID, INSIDE}},
{{12, 18, 13}, {INSIDE, INSIDE, INSIDE}},
{{13, 18, 19}, {INSIDE, RIGID, INSIDE}},
{{13, 19, 14}, {INSIDE, RIGID, INSIDE}},

{{20,  5,  6}, {RIGID, FREE, INSIDE}},
{{20,  6, 21}, {INSIDE, INSIDE, INSIDE}},
{{21,  6,  7}, {INSIDE, FREE, INSIDE}},
{{21,  7, 22}, {INSIDE, INSIDE, INSIDE}},
{{22,  7,  8}, {INSIDE, FREE, INSIDE}},
{{22,  8, 23}, {INSIDE, INSIDE, INSIDE}},
{{23,  8,  9}, {INSIDE, FREE, INSIDE}},
{{23,  9, 24}, {INSIDE, RIGID, INSIDE}},

{{25, 20, 21}, {RIGID, INSIDE, INSIDE}},
{{25, 21, 26}, {INSIDE, INSIDE, RIGID}},
{{26, 21, 22}, {INSIDE, INSIDE, INSIDE}},
{{26, 22, 27}, {INSIDE, INSIDE, RIGID}},
{{27, 22, 23}, {INSIDE, INSIDE, INSIDE}},
{{27, 23, 28}, {INSIDE, INSIDE, RIGID}},
{{28, 23, 24}, {INSIDE, INSIDE, INSIDE}},
{{28, 24, 29}, {INSIDE, RIGID, RIGID}}
};

Element_Data element_data[Element_Num];