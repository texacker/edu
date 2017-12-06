void init_data(void)
{
	int i, j, k;
	double swap;

	Element_Type *ep;
	Element_Data *ed;

	for(i = 0; i < Element_Num; i++)
	{
		ep = &element[i];
		ed = &element_data[i];

		for(j = 0; j < 3; j++)	//Get x, y of all three point.
		{
			ed->x[j] = ed->sx[j] = node[ep->n[j]].x;
			ed->y[j] = ed->sy[j] = node[ep->n[j]].y;
			ed->nx[j] = node[ep->n[(j+1)%3]].x;
			ed->ny[j] = node[ep->n[(j+1)%3]].y;
			ed->dx[j] = ed->nx[j] - ed->x[j];
			ed->dy[j] = ed->ny[j] - ed->y[j];
			ed->len[j] = sqrt(ed->dx[j] * ed->dx[j] + ed->dy[j] * ed->dy[j]);
			ed->nr[j] = ed->dy[j] / ed->len[j];
		}

		ed->xc = (ed->x[0] + ed->x[1] + ed->x[2]) / 3.0;
		ed->yc = (ed->y[0] + ed->y[1] + ed->y[2]) / 3.0;

		ed->a[0] = ed->x[1] * ed->y[2] - ed->x[2] * ed->y[1];
		ed->a[1] = ed->x[2] * ed->y[0] - ed->x[0] * ed->y[2];
		ed->a[2] = ed->x[0] * ed->y[1] - ed->x[1] * ed->y[0];

		ed->b[0] = ed->y[1] - ed->y[2];
		ed->b[1] = ed->y[2] - ed->y[0];
		ed->b[2] = ed->y[0] - ed->y[1];

		ed->c[0] = ed->x[2] - ed->x[1];
		ed->c[1] = ed->x[0] - ed->x[2];
		ed->c[2] = ed->x[1] - ed->x[0];

		ed->s = 0.5 * (ed->a[0] + ed->a[1] + ed->a[2]);

		for(j = 0; j < 3; j++)	//Sorting the point on x.
		{
			for(k = 0; k < 2; k++)
			{
				if (ed->sx[k] > ed->sx[k+1])
				{
					swap = ed->sx[k];
					ed->sx[k] = ed->sx[k+1];
					ed->sx[k+1] = swap;

					swap = ed->sy[k];
					ed->sy[k] = ed->sy[k+1];
					ed->sy[k+1] = swap;
		}}}

	}
}

void reset_Ae(void) { for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) Ae[i][j] = 0; }
void reset_Be(void) { for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) Be[i][j] = 0; }
void reset_Ce(void) { for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) Ce[i][j] = 0; }
void reset_De(void) { for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) De[i][j] = 0; }

void reset_L(void)
{
	for(int i = 0; i < Node_Num; i++)
		for(int j = 0; j < Node_Num; j++)
			A[i][j] = B[i][j] = C[i][j] = D[i][j] = 0;
}

void reset_G(void)
{
	for(int i = 0; i < 2 * (Node_Num + Free_Num); i++)
		for(int j = 0; j < 2 * (Node_Num + Free_Num); j++)
			G[i][j] = 0;
}

void gen_ACe(unsigned int e)
{
	int i, j, k, m;
	double eA, eB, eC;
	double t, r[2][2], z1[2][2], z2[2][2];

	Element_Data &ed = element_data[e];

	t = Gauss_T;
	r[0][0] = (t * (ed.sx[1] - ed.sx[0]) + (ed.sx[0] + ed.sx[1])) / 2;
	z1[0][0] = (ed.sy[2] - ed.sy[0]) * (ed.sx[1] - ed.sx[0]) * (t + 1)  / (ed.sx[2] - ed.sx[0]) / 2 + ed.sy[0];
	z2[0][0] = (t * (ed.sy[1] - ed.sy[0]) + (ed.sy[0] + ed.sy[1])) / 2;

	t = -Gauss_T;
	r[0][1] = (t * (ed.sx[1] - ed.sx[0]) + (ed.sx[0] + ed.sx[1])) / 2;
	z1[0][1] = (ed.sy[2] - ed.sy[0]) * (ed.sx[1] - ed.sx[0]) * (t + 1)  / (ed.sx[2] - ed.sx[0]) / 2 + ed.sy[0];
	z2[0][1] = (t * (ed.sy[1] - ed.sy[0]) + (ed.sy[0] + ed.sy[1])) / 2;

	t = Gauss_T;
	r[1][0] = (t * (ed.sx[1] - ed.sx[2]) + (ed.sx[2] + ed.sx[1])) / 2;
	z1[1][0] = (ed.sy[2] - ed.sy[0]) * (ed.sx[1] - ed.sx[2]) * (t + 1)  / (ed.sx[2] - ed.sx[0]) / 2 + ed.sy[2];
	z2[1][0] = (t * (ed.sy[1] - ed.sy[2]) + (ed.sy[2] + ed.sy[1])) / 2;

	t = -Gauss_T;
	r[1][1] = (t * (ed.sx[1] - ed.sx[2]) + (ed.sx[2] + ed.sx[1])) / 2;
	z1[1][1] = (ed.sy[2] - ed.sy[0]) * (ed.sx[1] - ed.sx[2]) * (t + 1)  / (ed.sx[2] - ed.sx[0]) / 2 + ed.sy[2];
	z2[1][1] = (t * (ed.sy[1] - ed.sy[2]) + (ed.sy[2] + ed.sy[1])) / 2;

	reset_Ae();
	reset_Ce();

	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		{
			Ce[i][j] = ed.c[i] * ed.c[j] * ed.xc / 4.0 / ed.s;
			Ae[i][j] = ed.b[i] * ed.b[j] * ed.xc / 4.0 / ed.s + Ce[i][j];

			t = 0;

			for(k = 0; k < 2; k++)
			{
				for(m = 0; m < 2; m++)
				{
					eA = (ed.a[i] + ed.b[i] * r[k][m]) * (ed.a[j] + ed.b[j] * r[k][m]);
					eB = ((ed.a[i] + ed.b[i] * r[k][m]) * ed.c[j] + (ed.a[j] + ed.b[j] * r[k][m]) * ed.c[i]) / 2;
					eC = (ed.c[i] * ed.c[j]) / 3;

					t += (
					eA * fabs(z2[k][m] - z1[k][m]) +
					eB * fabs(z2[k][m] * z2[k][m] - z1[k][m] * z1[k][m]) +
					eC * fabs(z2[k][m] * z2[k][m] * z2[k][m] - z1[k][m] * z1[k][m] * z1[k][m])) *
					(ed.sx[1] - ed.sx[k * 2]) / r[k][m];
			}}

			t *= Global_m * Global_m / 8 / ed.s / ed.s;

			Ae[i][j] += t;

	}}
}

void gen_Be(unsigned int e)
{
	int i, j, k;
	double t, r, z;
	double tmp;

	Element_Type &ep = element[e];
	Element_Data &ed = element_data[e];

	reset_Be();

	for(i = 0; i < 3; i++)
	{
		if ((ep.e[i] == FREE)||(ep.e[i] == RIGID))
		{
			for(j = 0; j < 3; j++)
				for(k = 0; k < 3; k++)
				{
					t = Gauss_T;
					r = (t * (ed.nx[i] - ed.x[i]) + (ed.nx[i] + ed.x[i])) / 2;
					z = (t * (ed.ny[i] - ed.y[i]) + (ed.ny[i] + ed.y[i])) / 2;
					tmp = (ed.a[j] + ed.b[j] * r + ed.c[j] * z) * (ed.a[k] + ed.b[k] * r + ed.c[k] * z);

					t = -Gauss_T;
					r = (t * (ed.nx[i] - ed.x[i]) + (ed.nx[i] + ed.x[i])) / 2;
					z = (t * (ed.ny[i] - ed.y[i]) + (ed.ny[i] + ed.y[i])) / 2;
					tmp += (ed.a[j] + ed.b[j] * r + ed.c[j] * z) * (ed.a[k] + ed.b[k] * r + ed.c[k] * z);

					Be[j][k] += tmp * ed.nr[i] * ed.len[i] / 8 / ed.s / ed.s;
	}}}
}

void gen_De(unsigned int e)
{
	int i, j, k;
	double t, r, z;
	double tmp;

	Element_Type &ep = element[e];
	Element_Data &ed = element_data[e];

	reset_De();

	for(i = 0; i < 3; i++)
	{
		if (ep.e[i] == FREE)
		{
			for(j = 0; j < 3; j++)
				for(k = 0; k < 3; k++)
				{
					t = Gauss_T;
					r = (t * (ed.nx[i] - ed.x[i]) + (ed.nx[i] + ed.x[i])) / 2;
					z = (t * (ed.ny[i] - ed.y[i]) + (ed.ny[i] + ed.y[i])) / 2;
					tmp = (ed.a[j] + ed.b[j] * r + ed.c[j] * z) * (ed.a[k] + ed.b[k] * r + ed.c[k] * z) / sqrt(1 + (Fea_Len / r) * (Fea_Len / r));

					t = -Gauss_T;
					r = (t * (ed.nx[i] - ed.x[i]) + (ed.nx[i] + ed.x[i])) / 2;
					z = (t * (ed.ny[i] - ed.y[i]) + (ed.ny[i] + ed.y[i])) / 2;
					tmp += (ed.a[j] + ed.b[j] * r + ed.c[j] * z) * (ed.a[k] + ed.b[k] * r + ed.c[k] * z) / sqrt(1 + (Fea_Len / r) * (Fea_Len / r));

					De[j][k] += tmp * ed.nr[i] * ed.len[i] / 8 / ed.s / ed.s;
	}}}
}

void gen_L(void)
{
	int i, j, k;
	Element_Type *ep;
	Element_Data *ed;

	reset_L();

	for(i = 0; i < Element_Num; i++)
	{
		ep = &element[i];
		ed = &element_data[i];

		gen_ACe(i);
		gen_Be(i);
		gen_De(i);

		for(j = 0; j < 3; j++)
		{
			for(k = 0; k < 3; k++)
			{
				A[ep->n[j]][ep->n[k]] += Ae[j][k] + De[j][k] * 4;
				B[ep->n[j]][ep->n[k]] += Be[j][k];
				C[ep->n[j]][ep->n[k]] += Ce[j][k];
				D[ep->n[j]][ep->n[k]] += De[j][k] * 4;
			}
		}
	}
}

void gen_total(void)
{
	int i, j, u, v;
	int l = Node_Num;
	int m = Free_Num;
	int n = l - m;

	double *d11 = (double *)malloc((m * m) * sizeof(double));
	double *d11I = (double *)malloc((m * m) * sizeof(double));

	double *a11 = (double *)malloc((m * m) * sizeof(double));
	double *a12 = (double *)malloc((m * n) * sizeof(double));
	double *a21 = (double *)malloc((m * n) * sizeof(double));
	double *a22 = (double *)malloc((n * n) * sizeof(double));
	double *a11E = (double *)malloc((m * m) * sizeof(double));
	double *a21S = (double *)malloc((m * n) * sizeof(double));
	double *a22I = (double *)malloc((n * n) * sizeof(double));

	double *b11 = (double *)malloc((m * m) * sizeof(double));
	double *b12 = (double *)malloc((m * n) * sizeof(double));
	double *b21 = (double *)malloc((m * n) * sizeof(double));
	double *b22 = (double *)malloc((n * n) * sizeof(double));
	double *b11E = (double *)malloc((m * m) * sizeof(double));
	double *b12E = (double *)malloc((m * n) * sizeof(double));
	double *b21S = (double *)malloc((m * n) * sizeof(double));
	double *b22S = (double *)malloc((n * n) * sizeof(double));

	double *c11 = (double *)malloc((m * m) * sizeof(double));
	double *c12 = (double *)malloc((m * n) * sizeof(double));
	double *c21 = (double *)malloc((m * n) * sizeof(double));
	double *c22 = (double *)malloc((n * n) * sizeof(double));
	double *c11E = (double *)malloc((m * m) * sizeof(double));
	double *c12E = (double *)malloc((m * n) * sizeof(double));
	double *c21S = (double *)malloc((m * n) * sizeof(double));
	double *c22S = (double *)malloc((n * n) * sizeof(double));

	reset_G();
	gen_L();

	for (i = 0; i < (Node_Num + Free_Num); i++) G[i][Node_Num+i] = 1;
	for (i = 0; i < Free_Num; i++) G[2 * Node_Num + i][2 * Node_Num + Free_Num + i] = 1;

	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
		{
			u = i * m + j;
			a11[u] = A[i][j];
			b11[u] = B[i][j];
			c11[u] = C[i][j];
			d11[u] = d11I[u] = D[i][j];
		}

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			u = i * n + j;
			v = j * m + i;
			a12[u] = A[i][j+m];
			b12[u] = B[i][j+m];
			c12[u] = C[i][j+m];

			a21[v] = A[j+m][i];
			b21[v] = B[j+m][i];
			c21[v] = C[j+m][i];
		}

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			u = i * n + j;
			a22[u] = a22I[u] = A[i+m][j+m];
			b22[u] = B[i+m][j+m];
			c22[u] = C[i+m][j+m];
		}

	rinv(a22I, n);
	rinv(d11I, m);

	rmul(a12, rmul(a22I, a21, n, n, m, a21S), m, n, m, a11E);
	rmul(a12, rmul(a22I, b21, n, n, m, b21S), m, n, m, b11E);
	rmul(a12, rmul(a22I, b22, n, n, n, b22S), m, n, n, b12E);
	rmul(a12, rmul(a22I, c21, n, n, m, c21S), m, n, m, c11E);
	rmul(a12, rmul(a22I, c22, n, n, n, c22S), m, n, n, c12E);

	radd(a11, rmul(a11E, a11E, m, m, -1), m, m, a11E);
	radd(b11, rmul(b11E, b11E, m, m, -1), m, m, b11E);
	radd(b12, rmul(b12E, b12E, m, n, -1), m, n, b12E);
	radd(c11, rmul(c11E, c11E, m, m, -1), m, m, c11E);
	radd(c12, rmul(c12E, c12E, m, n, -1), m, n, c12E);

	rmul(d11I, a11E, m, m, m, a11E);
	rmul(d11I, b11E, m, m, m, b11E);
	rmul(d11I, b12E, m, m, n, b12E);
	rmul(d11I, c11E, m, m, m, c11E);
	rmul(d11I, c12E, m, m, n, c12E);

	for (i = 0; i < m; i++)
	{
		v = 2 * l + m + i;
		for (j = 0; j < m; j++)
		{
			u = i * m + j;
			G[v][j] = c11E[u] / -4;
			G[v][l+j] = b11E[u] * Global_m / 4;
			G[v][2*l+j] = a11E[u] / 4;
		}
		G[v][2*l+i] += 1;
	}

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			u = i * n + j;
			v = j * m + i;
			G[2*l+m+i][m+j] = -c12E[u] / 4;
			G[2*l+m+i][l+m+j] = b12E[u] * Global_m / 4;
			G[l+m+j][i] = c21S[v];
			G[l+m+j][l+i] = -b21S[v] * Global_m;
			G[l+m+j][2*l+i] = -a21S[v];
		}

	for (i = 0; i < n; i++)
	{
		v = l + m + i;
		for (j = 0; j < n; j++)
		{
			u = i * n + j;
			G[v][m+j] = c22S[u];
			G[v][l+m+j] = -b22S[u] * Global_m;
		}
	}

	free(d11);
	free(d11I);

	free(a11);
	free(a12);
	free(a21);
	free(a22);
	free(a11E);
	free(a21S);
	free(a22I);

	free(b11);
	free(b12);
	free(b21);
	free(b22);
	free(b11E);
	free(b12E);
	free(b21S);
	free(b22S);

	free(c11);
	free(c12);
	free(c21);
	free(c22);
	free(c11E);
	free(c12E);
	free(c21S);
	free(c22S);
}