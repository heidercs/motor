#ifndef BANDDTW_H
#define BANDDTW_H

#include "lineal.h"

enum TYPE_BAND_DTW { NON_BAND=0, SAKOE_CHIBA=1, ITAKURA=2} ;

/** support multidimensional objects  **/

class BandDTW: public Distance
{
private:    
    // it belongs to TYPE_BAND_DTW
    int type_band;
    
    // fixed the maximum radio of the band
    int r_max;    

    // Also, the r_max can be update in runtime 
    // applying a fraction of the total size  n  
    // Band = 2*r_max && r_max = fraction_r * n
    float fraction_r;

    // select the base distance function
    Lineal* dist_lineal;
    

public:

    BandDTW(int _type_band=NON_BAND, double _fraction_r=0.25, int _type_lineal=2);
    void set_band(int _type_band, double _fraction_r);
    void set_band(int _type_band, int _r_max);

    int get_type_lineal();
    
	double* get_U(double* q,int nq);
    double* get_L(double* q,int nq);        

    double d(Object* p, Object* q);
    double lower_bound(Object* p, Object* q);
    double upper_bound(Object* p, Object* q);
    
    double lower_bound(Object *p, double *L, double *U, int n);
    double upper_bound(Object *p, double *L, double *U, int n);

    // we use anticipatory pruning for computing the distance by levels
    double d(Object* p, Object* q, double th_dist); 

private:
    int bandMin(int pos, int dimension);
    int bandMax(int pos, int dimension);
    void fix_r_max(int n);
    double dist(double *p, double *q, int np, int nq);
    double dist(double **p, double **q,  int np, int nq, int D);
    double dist(double **p, double **q, double *wp, double *wq, int np, int nq, int D);
    double dist_by_levels(double **p, double **q,  int np, int nq, int D, double th_dist = DBL_MAX);
    double partial_dist(double **p, double **q,  int izq, int der, int D);

    double dist(double *vp, double *vq, double *sp, double *sq, int np, int nq);    
    double dist_by_levels(double* vp, double* vq, double* sp, double* sq,  int np, int nq, double th_dist = DBL_MAX);       
    double partial_dist(double *vp, double *vq, double *sp, double *sq, int izq, int der); 

    double dist(double **vp, double **vq, double **sp, double **sq, int np, int nq, int D);    
    double dist_by_levels(double** vp, double** vq, double** sp, double** sq,  int np, int nq, int D, double th_dist = DBL_MAX);       
    double partial_dist(double** vp, double** vq, double** sp, double** sq, int izq, int der, int D); 

    
   
};

#endif // BANDADTW_H
