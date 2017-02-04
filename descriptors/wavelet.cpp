#include "wavelet.h"

double* Wavelet::get_haar(double *vec, int n, int new_n)
{
    int i=0;
    double *vecp = new double[n];
    double *haar = new double[n];
    int w=n;

    for(i=0; i<n; i++)
    {
        vecp[i] = 0;
        haar[i]=vec[i];
    }

    double by2 = 2.0; // originalmente se usa sqrt(2.0)

    while(w>1)
    {
        w/=2;
        for(i=0; i<w; i++)
        {
            vecp[i] = (haar[2*i] + haar[2*i+1])/by2;
            vecp[i+w] = (haar[2*i] - haar[2*i+1])/by2;
        }

        for(i=0; i<(w*2); i++)
            haar[i] = vecp[i];
    }

    delete[] vecp;

    if(new_n < n) // keep some coefficients
	{
	    double *new_vec = new double[new_n];
	    for(i=0;i<new_n;i++)
	    	new_vec[i]=haar[i];

	    delete[] haar;

	    return new_vec;
	}

    return haar;
}


/** The 2D Haar Transform **/
double** Wavelet::get_haar2D(double **matrix, int n, int m)
{
	double *temp_row = new double[m];
	double *temp_col = new double[n];
	double **haar = new double*[n];
	double *temp;

	int i=0,j=0;
	int w = m, h=n;

	for(i=0; i<n; i++)
    {
        haar[i]=new double[m];
        for(j=0; j<m; j++)
            haar[i][j]=matrix[i][j];
    }

	while(w>1 || h>1)
	{
		if(w>1)
		{
			for(i=0;i<h;i++)
			{
				for(j=0;j<m;j++)
					temp_row[j] = haar[i][j];

				temp=get_haar(temp_row,m,w);
				delete[] temp_row;
				temp_row=temp;

				for(j=0;j<m;j++)
					haar[i][j] = temp_row[j];
			}
		}

		if(h>1)
		{
			for(i=0;i<w;i++)
			{
				for(j=0;j<n;j++)
					temp_col[j] = haar[j][i];
				
				temp=get_haar(temp_col, n, h);
				delete[] temp_col;
				temp_col=temp;
				
				for(j=0;j<n;j++)
					haar[j][i] = temp_col[j];
			}
		}

		if(w>1)
			w/=2;
		if(h>1)
			h/=2;
	}

	delete[] temp_row;
	delete[] temp_col;

	return haar;
}
