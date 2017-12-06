//This is the beganning of my thesis which started at May 25, 1998

#include <iostream.h>
#include <conio.h>
#include <math.h>
#include <alloc.h>
#include <stdio.h>
#include "struct.h"
#include "math_lib.h"
#include "thesis_1.h"	//#include "thesis_0.h"

void main(void)
{
	int i, j, k, jt = 60;
	Element_Type *ep;
	Element_Data *ed;

	clrscr();
	cout << "Call init porcessure." << endl;

	init_data();
	gen_total();
	chhbg(&G[0][0], 2 * (Node_Num + Free_Num));
	chhqr(&G[0][0], 2 * (Node_Num + Free_Num), vr, vi, Eps, jt);

	for(i = 0; i < 8; i++)
	{
//		for(j = 0; j < 2 * (Node_Num + Free_Num); j++)
//		{
			cout << vr[i] << " " << vi[i] << " ";
//		}
		cout << endl;
	}

/*	for(i = 0; i < Node_Num; i++)
	{
		for(j = 0; j < Node_Num; j++)
		{
			cout << C[i][j] << " ";
		}
		cout << endl;
	}*/

//	cout << "The end of testing." << endl;
}