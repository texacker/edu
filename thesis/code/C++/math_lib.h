double *radd(double *a, double *b, int m, int k, double *c)
{
	int i, j, u;

	for (i = 0; i < m; i++)
		for (j = 0; j < k; j++)
		{
			u = i * k + j;
			c[u] = a[u] + b[u];
		}
	return c;
}

double *rmul(double *a, double *b, int m, int k, double c)
{
	int i, j, u;

	for (i = 0; i < m; i++)
		for (j = 0; j < k; j++)
		{
			u = i * k + j;
			b[u] = a[u] * c;
		}
	return b;
}

double *rmul(double *a, double *b, int m, int n, int k, double *c)
{
	int i, j, l, u;
	double *d = (double *)malloc((m * k) * sizeof(double));

	for (i = 0; i < m; i++)
		for (j = 0; j < k; j++)
		{
			u = i * k + j;
			d[u] = 0.0;
			for (l = 0; l < n; l++) d[u] += a[i * n + l] * b[l * k + j];
			c[u] = d[u];
		}
	free(d);
	return c;
}

int rinv(double *a, int n)
{
	int i, j, k, l, u, v;
	double d, p;
	int *is = (int *)malloc(n * sizeof(int));
	int *js = (int *)malloc(n * sizeof(int));

	for (k = 0; k <= n-1; k++)
	{
		d = 0.0;
		for (i = k; i <= n-1; i++)
			for (j = k; j <= n-1; j++)
			{
				l =i * n + j;
				p = fabs(a[l]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}

		if ((d + 1.0) == 1.0)
		{
			free(is);
			free(js);
			cerr << "err**not inv" << endl;
			return(0);
		}

		if (is[k] != k)
			for (j = 0; j <= n-1; j++)
			{
				u = k * n + j;
				v = is[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}

		if (js[k] != k)
			for (i = 0; i <= n-1; i++)
			{
				u =i * n + k;
				v = i * n + js[k];
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}

		l = k * n + k;
		a[l] = 1.0 / a[l];

		for (j = 0; j <= n-1; j++)
			if (j != k)
			{
				u = k * n + j;
				a[u] = a[u] * a[l];
			}
		for (i = 0; i <= n-1; i++)
			if (i != k)
				for (j = 0; j <= n-1; j++)
					if (j != k)
					{
						u = i * n + j;
						a[u] = a[u] - a[i*n+k] * a[k*n+j];
					}

		for (i = 0; i <= n-1; i++)
			if (i != k)
			{
				u = i * n + k;
				a[u] = -a[u] * a[l];
			}
	}

	for (k = n-1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n-1; j++)
			{
				u = k * n + j;
				v = js[k] * n + j;
				p = a[u];
				a[u] = a[v];
				a[v] = p;
			}

		if (is[k] != k)
			for (i = 0; i <= n-1; i++)
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
	return(1);
}

void chhbg(double *a, int n)
{
	int i, j, k, u, v;
	double d, t;

	for (k = 1; k <= n-2; k++)
	{
		d = 0.0;
		for (j = k; j <= n-1; j++)
		{
			u = j * n + k-1;
			t = a[u];
			if (fabs(t) > fabs(d))
			{
				d = t;
				i = j;
			}
		}
		if (fabs(d)+1.0 != 1.0)
		{
			if (i != k)
			{
				for (j = k-1; j <= n-1; j++)
				{
					u = i * n + j;
					v = k * n + j;
					t = a[u];
					a[u] = a[v];
					a[v] = t;
				}
				for (j = 0; j <= n-1; j++)
				{
					u = j * n + i;
					v = j * n + k;
					t = a[u];
					a[u] = a[v];
					a[v] = t;
				}
			}
			for (i = k+1; i <= n-1; i++)
			{
				u = i * n + k-1;
				t = a[u] / d;
				a[u] = 0.0;
				for (j = k; j <= n-1; j++)
				{
					v = i * n + j;
					a[v] = a[v] - t * a[k*n+j];
				}
				for (j = 0; j <= n-1; j++)
				{
					v = j * n + k;
					a[v] = a[v] + t * a[j*n+i];
				}
			}
		}
	}
	return;
}

int chhqr(double *a, int n, double *u, double *v, double eps, int jt)
{
	int m, it, i, j, k, l, ii, jj, kk, ll;
	double b, c, w, g, xy, p, q, r, x, s, e, f, z, y;
	it = 0;
	m = n;
	while (m != 0)
	{
		l = m - 1;
		while ((l>0) && (fabs(a[l*n+l-1]) > eps*(fabs(a[(l-1)*n+l-1])+fabs(a[l*n+l])))) l = l - 1;
		ii = (m - 1) * n + m - 1;
		jj = (m - 1) * n + m - 2;
		kk = (m - 2) * n + m - 1;
		ll = (m - 2) * n + m - 2;
		if (l == m-1)
		{
			u[m-1] = a[(m-1)*n+m-1];
			v[m-1] = 0.0;
			m = m-1;
			it = 0;
		}
		else
		if (l == m-2)
		{
			b = -(a[ii] + a[ll]);
			c = a[ii] * a[ll] - a[jj] * a[kk];
			w = b * b - 4.0 * c;
			y = sqrt(fabs(w));
			if (w > 0.0)
			{
				xy = 1.0;
				if (b < 0.0) xy = -1.0;
				u[m-1] = (-b - xy * y) / 2.0;
				u[m-2] = c / u[m-1];
				v[m-1] = 0.0;
				v[m-2] = 0.0;
			}
			else
			{
				u[m-1] = -b / 2.0;
				u[m-2] = u[m-1];
				v[m-1] = y / 2.0;
				v[m-2] = -v[m-1];
			}
			m = m - 2;
			it = 0;
		}
		else
		{
			if (it >= jt)
			{
				printf("fail\n");
				return(-1);
			}
			it = it + 1;
			for (j = l+2; j <= m-1; j++)
				a[j*n+j-2] = 0.0;
			for (j = l+3; j <= m-1; j++)
				a[j*n+j-3] = 0.0;
			for (k = l; k <= m-2; k++)
			{
				if (k != l)
				{
					p = a[k*n+k-1];
					q = a[(k+1)*n+k-1];
					r = 0.0;
					if (k != m-2) r = a[(k+2)*n+k-1];
				}
				else
				{
					x = a[ii] + a[ll];
					y = a[ll] * a[ii] - a[kk] * a[jj];
					ii = l * n + l;
					jj = l * n + l + 1;
					kk = (l+1)*n+l;
					ll = (l+1)*n+l+1;
					p = a[ii] * (a[ii] - x) + a[jj] * a[kk] + y;
					q = a[kk] * (a[ii] + a[ll] - x);
					r = a[kk] * a[(l+2)*n+l+1];
				}
				if ((fabs(p) + fabs(q) + fabs(r)) != 0.0)
				{
					xy = 1.0;
					if (p < 0.0) xy = -1.0;
					s = xy * sqrt(p*p+q*q+r*r);
					if (k != l) a[k*n+k-1] = -s;
					e = -q / s;
					f = -r / s;
					x = -p / s;
					y = -x - f * r / (p + s);
					g = e * r / (p + s);
					z = -x - e * q / (p + s);
					for (j = k; j <= m-1; j++)
					{
						ii = k * n + j;
						jj= (k + 1) * n + j;
						p = x * a[ii] + e * a[jj];
						q = e * a[ii] + y * a[jj];
						r = f * a[ii] + g * a[jj];
						if (k != m-2)
						{
							kk = (k+2) * n + j;
							p = p + f * a[kk];
							q = q + g * a[kk];
							r = r + z * a[kk];
							a[kk] = r;
						}
						a[jj] = q;
						a[ii] = p;
					}
					j = k + 3;
					if (j >= m-1) j = m-1;
					for (i = l; i <= j; i++)
					{
						ii = i * n + k;
						jj = i * n + k + 1;
						p = x * a[ii] + e * a[jj];
						q = e * a[ii] + y * a[jj];
						r = f * a[ii] + g * a[jj];
						if (k != m-2)
						{
							kk = i * n + k + 2;
							p = p + f * a[kk];
							q = q + g * a[kk];
							r = r + z * a[kk];
							a[kk] = r;
						}
						a[jj] = q;
						a[ii] = p;
					}
				}
			}
		}
	}
	return(1);
}