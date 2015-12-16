#include "../objects/vector.h"
#include "../objects/matrix.h"
#include "../objects/vectorpair.h"
#include "../objects/matrixpair.h"
#include "../descriptors/sax.h"
#include "../util/functions.h"

///////////////////////////////////////////////////////////////////////////////////////////
 /**
  * [Apply the TVA-representation to the time series]
  * @param  p       [data vector]
  * @param  N       [vector length]
  * @param  n       [total number of segments]
  */ 
inline VectorPair* get_vectorpair(Vector* P, int n)
{
	double* p = P->get_elems();
	int N = P->get_n();
	n = std::min(N/2, n); // adjust the user segments number 

	double w = 1.0*N/n;
	double* coef;
	int init, end;

	double *slope = new double[n];
	double *value = new double[n];
	double *x = new double[(int)ceil(w)];
	double *y = new double[(int)ceil(w)];	
	    
	init = 0;
	for (int i = 0; i < n; ++i)
	{
		end=std::min((int) round((i+1) * w), N);
	
		for (int j = init; j < end; ++j)
		{
			x[j-init] = (j + 1) * 2.0 / w; 
			y[j-init] = p[j]; 
		}

		coef =  lm(x, y, end - init);

		slope[i] = coef[0]; // slope of the segment
		//slope[i] = coef[0] / 1.5707963267948966; // normalized slope [-1 : +1]
		value[i] = coef[2]; // mean value		

		delete[] coef;

		init = end;
	}

	delete[] x;
	delete[] y;


	return new VectorPair(value, slope, n, P->get_id_file());
}


/**
 * [Apply the TVA-representation to the time series hierarchical]
 * @param  p       [data vector]
 * @param  N       [vector length]
 * @param  L       [resolution levels > 0]
 */
inline VectorPair* get_hie_vectorpair(Vector* P, int L = -1, bool by_levels=false)
{
	double* p = P->get_elems();
	int N =  P->get_n();

	//int Lmax = log2(N/2.0  +  1); // max level of p
	int Lmax = log2(N/2.0); // max level of p
	if(L > 0)
		L = std::min(Lmax, L); // adjust the user level
	else
		L = Lmax;
	
	int total = pow(2, L + 1) - 1; // total number of pairs

	double *x, *y, *coef, w;
	int n, init, end, i,j,l, k = 0;

	double *slope = new double[total];
	double *value = new double[total];
		
	for (l = 0; l <= L; ++l)
	{
		n = pow(2, l);
		w = 1.0*N/n;

		x = new double[(int)ceil(w)];
		y = new double[(int)ceil(w)];

		init = 0;
		for (i = 0; i < n; ++i)
		{
			end=std::min((int) round((i+1) * w), N);
		
			for (j = init; j < end; ++j)
			{
				x[j-init] = (j + 1)* 2.0 / w; // normalize in each level
				y[j-init] = p[j]; 
			}

			coef =  lm(x, y, end - init);

			slope[k] = coef[0]; // slope of the segment
			//slope[k] = coef[0] / 1.5707963267948966; // normalized slope [-1 : +1]
			value[k] = coef[2]; // mean value

			delete[] coef;

			init = end;
			k++;
		}

		delete[] x;
		delete[] y;
	}

	return new VectorPair(value, slope, total, P->get_id_file(), by_levels);	
}


/**
 * [sax_vectorpair description]
 * @param  P         [description]
 * @param  value_sax [description]
 * @param  slope_sax [description]
 * @return           [description]
 */
inline VectorPair* sax_vectorpair(VectorPair* P, Sax *value_sax, Sax* slope_sax)
{
	int n = P->get_n();
	int *v_sax = value_sax->run_sax_from_paa(P->get_value(), n);
    int *s_sax = slope_sax->run_sax_from_paa(P->get_slope(), n);

    double *v_double = new double[n];
    double *s_double = new double[n];
    std::copy(v_sax, v_sax + n, v_double);
    std::copy(s_sax, s_sax + n, s_double);

    delete[] v_sax;
    delete[] s_sax;

    return new VectorPair(v_double, s_double, n, P->get_id_file());
}


///////////////////////////////////////////////////////////////////////////////////////////
/**
* [Apply the TVA-representation to the time series]
* @param  Pinv    [data matrix]
* @param  n       [total number of segments]
*/

inline MatrixPair* get_matrixpair(Matrix* Pinv, int n)
{
	double** ma = Pinv->get_elems();
	int D = Pinv->get_n(); // variables
	int N = Pinv->get_D(); // matrix length
	n = std::min(N/2, n); // adjust the user segments number 

	double w = 1.0*N/n;
	double* coef;
	int init, end, i, j, d;

	double **slope = new double*[D];
	double **value = new double*[D];	
	double *x = new double[(int)ceil(w)];
	double *y = new double[(int)ceil(w)];	

	for (d = 0; d < D; d++)
    {
        slope[d] = new double[n];
        value[d] = new double[n];
    }
	    
	init = 0;
	for (i = 0; i < n; ++i)
	{
		end=std::min((int) round((i+1) * w), N);
	
		for (d = 0; d < D; ++d)
		{			
		
			for (j = init; j < end; ++j)
			{
				x[j-init] = (j + 1) * 2.0 / w; 
				y[j-init] = ma[d][j]; 
			}

			coef =  lm(x, y, end - init);

			slope[d][i] = coef[0]; // slope of the segment
			//slope[i] = coef[0] / 1.5707963267948966; // normalized slope [-1 : +1]
			value[d][i] = coef[2]; // mean value		

			delete[] coef;
		}

		init = end;
	}

	double **value_t = transposeT(value, D, n);
	double **slope_t = transposeT(slope, D, n);

	///////////// free memory ////////////////
	for (int i = 0; i < D; ++i)
	{
		delete[] value[i];	
		delete[] slope[i];
	}
	delete[] value;	
	delete[] slope;
	delete[] x;
	delete[] y;
	////////////////////////////////////////

	return new MatrixPair(value_t, slope_t, n, D, Pinv->get_id_file());
}




/**
 * [Apply the TVA-representation to the multivariate time series hierarchical]
 * @param  Pinv     [data matrix]
 * @param  L 	    [resolution levels > 0]
 */
inline MatrixPair* get_hie_matrixpair(Matrix* Pinv, int L = -1, bool by_levels=false)
{
	double** ma = Pinv->get_elems();
	int D = Pinv->get_n(); // variables
	int N = Pinv->get_D(); // matrix length

	int Lmax = log2(N/2.0); // max level of p
	if(L > 0)
		L = std::min(Lmax, L); // adjust the user level
	else
		L = Lmax;
	
	int total = pow(2, L + 1) - 1; // total number of pairs

	double *x, *y, *coef, w;
	int n, init, end, i,j,l, d, k = 0;

	double **slope = new double*[D];
	double **value = new double*[D];

	for (d = 0; d < D; d++)
    {
        slope[d] = new double[total];
        value[d] = new double[total];
    }
		
	for (l = 0; l <= L; ++l)
	{
		n = pow(2, l);
		w = 1.0*N/n;

		x = new double[(int)ceil(w)];
		y = new double[(int)ceil(w)];

		init = 0;
		for (i = 0; i < n; ++i)
		{
			end=std::min((int) round((i+1) * w), N);

			for (d = 0; d < D; ++d)
			{
		
				for (j = init; j < end; ++j)
				{
					x[j-init] = (j + 1)* 2.0 / w; // normalize in each level
					y[j-init] = ma[d][j]; 
				}

				coef =  lm(x, y, end - init);

				slope[d][k] = coef[0]; // slope of the segment
				//slope[k] = coef[0] / 1.5707963267948966; // normalized slope [-1 : +1]
				value[d][k] = coef[2]; // mean value					

				delete[] coef;
			}

			init = end;
			k++;
		}

		delete[] x;
		delete[] y;
	}


	double **value_t = transposeT(value, D, total);
	double **slope_t = transposeT(slope, D, total);

	///////////// free memory ////////////////
	for (int i = 0; i < D; ++i)
	{
		delete[] value[i];	
		delete[] slope[i];
	}
	delete[] value;	
	delete[] slope;
	////////////////////////////////////////

	return new MatrixPair(value_t, slope_t, total, D, Pinv->get_id_file(), by_levels);	
}

///////////////////////////////////////////////////////////////////////////////
