#include "../descriptors/paa.h"
#include "../descriptors/pip.h"
#include "../descriptors/fourier.h"
#include "../descriptors/wavelet.h"
#include "../descriptors/sax.h"
#include "../objects/matrix.h"
#include "../util/functions.h"


/* return the name of descriptor  */
inline string get_name_descriptor(int type_descr)
{
    string name_descr;
    switch(type_descr)
    {
    case 0: name_descr = "stats"; break;   

    case 1: name_descr = "fourier"; break;
    
    case 2: name_descr = "wavelet"; break;

    case 3: name_descr = "paa"; break;

    case 4: name_descr = "pip"; break;    

    case 5: name_descr = "sax"; break; 

    default:  name_descr="none"; break;      
    }

    return name_descr;
}

/* reduce the dimensionality from n to new_n using any descriptor */
inline double* apply_descriptor(double *window, int N, int new_n, int type)
{
	double *coefs=NULL;
	
	switch(type)
    {    
    case 1: // Fourier
        coefs = Fourier::get_coef_fourier(window, N, new_n);
        break;
    case 2: // Haar Wavelet
        coefs = Wavelet::get_haar(window, N, new_n);
        break;
    case 3: // PAA
        coefs = Paa::get_paa(window, N, new_n);
        break;
    case 4: // PIP
        coefs = Pip::get_seg_pip(window, N, new_n);
        break;    
    case 5: // Sax
        coefs = Sax::get_sax_double(window, N, new_n);
        break;    
    default:
        std::cerr<<"This descriptor is not found\n";
        break;               
    }
	
	return coefs;
}

/* Apply any feature extraction technique to each component of the window  
 * Feature is conformed by spactral, time, singualarities or symbolic coefficients*/ 
inline double* apply_descriptor(double** window, int D,  int N, int new_n, int &new_dim, int type_descr = 1)
{       
    new_dim = D * new_n;
        
    double* vec = new double[new_dim];
    double* coefs = NULL;

    int i, j, k = 0;

    for (i = 0; i < D; ++i)
    {
        coefs=apply_descriptor(window[i], N, new_n, type_descr);

        for(j = 0; j < new_n; j++)
            vec[k++] = coefs[j];
        
        delete[] coefs;
     }

    return vec;
}

/* Feature is conformed by the Mean, Standar Deviation and Correlation between variables */
inline double* apply_stats(double** window, int D,  int N, int &new_dim)
{
    // mean(X,Y,Z), var(X,Y,Z), cor(X,Y), cor(Y,Z), cor(X, Z)
    new_dim = D + D + D*(D-1)/2; 

    double* vec = new double[new_dim];

    int k = 0;
    for (int i = 0; i < D; ++i)
    {   
        vec[k++] = media(window[i], N);                      
        vec[k++] = desviacion(window[i], N);                  
    }

    for (int i = 0; i < D - 1; ++i)
    {
        for (int j = i+1; j < D; ++j)
        {       
            vec[k++] = correlation(window[i], window[j], N);                  
        }
    }

    return vec;
}


//////////////// Global Feature - UTS//////////////////////

Vector* get_global(Vector* P, int &new_n, int type_descr)
{
    double* vec = P->get_elems();
    int N = P->get_n();    
    int id_file = P->get_id_file();   
           
    double *feature = apply_descriptor(vec, N, new_n, type_descr);    

    return new Vector(feature, new_n, id_file);        
}


//////////////// Global Feature - MTS//////////////////////
Object* get_global(Matrix* Pinv, int new_n, int type_descr)
{
    double** ma = Pinv->get_elems();
    int D = Pinv->get_n();
    int N = Pinv->get_D();
    int id_file = Pinv->get_id_file();

    Object* Q = NULL;

    if(type_descr == 0)                
    {
        double* coefs = apply_stats(ma, D, N, new_n);
        Q = new Vector(coefs, new_n, id_file);
    }
    else
    {
        double** coefs= new double*[D];
        
        new_n = std::min(new_n, N);

        for (int i = 0; i < D; ++i)
        {            
            coefs[i] = apply_descriptor(ma[i], N, new_n, type_descr);
        }

        Q = new Matrix(coefs, D, new_n, id_file, false, true);
        ((Matrix*)Q)->transpose();
    }

    return Q;
}