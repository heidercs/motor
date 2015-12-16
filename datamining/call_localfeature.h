#include <cmath>
#include <float.h>
#include <string.h>

#include "../objects/signature.h"
#include "../datamining/kmeans.h"
#include "../datamining/kmedoids.h"
#include "../distances/lineal.h"
#include "call_globalfeature.h"


////////////////////// Get Local Features from UTS ///////////////////////////

 /* 
 * @param  P             [data vector]
 * @param  nw            [size of the sliding window]
 * @param  overlap       [size of the overlap]
 * @param  n             [number total of feature vectors] 
 * @param  dim           [dimension of feature space]
 * @param  type_descr    [type_descr of the descriptor]
 * @param  save_position [if we want to save the referential position]
 * */  

/**
 * Get the signature receiving as parameter the sliding windows and the overlap,
 * The weight is equals to one
 */
Signature* get_signature(Vector* P, int nw, int overlap, int dim, int type_descr, bool save_position=false)
{
    int i, j, n, init, end, dim_add;    
    double **features, *window, *coefs;

    double* p = P->get_elems();
    int N = P->get_n(); //vector length
        
    // number of feature vectors
    n=ceil((N-nw)/(1.0*overlap))+1;
    // dim must not be greater than the size of window
    dim = std::min(dim, nw); 
    
    features=new double*[n];  
    
    for(i = 0; i < n; i++)
    {
        init = i*overlap;

        end = std::min(init + nw, N);
        window=new double[end - init];
        for (j = init; j < end  ; j++)      
            window[j - init] = p[j];
        
        features[i]=new double[dim];  
        
        dim_add =  save_position? 1 : 0;     
               
        coefs = apply_descriptor(window, end - init, dim - dim_add, type_descr);
                        
        for(j = dim_add; j < dim ; j++)
        {
            features[i][j] = coefs[j - dim_add];                
        }

        if(save_position)// the time series must be normalized, gaussian distribution
        {// adjust the position  between  [-1  :  1],       
            features[i][0] = equivalence(i, 0, n - 1, -1, 1);
        }       

        delete[] window;
        delete[] coefs;
    }

    return new Signature(features, NULL, n, dim, P->get_id_file());    
}

/** Get the signature receiving as parameter the number of feature vectors,
**  The weight is equals to one */
Signature* get_signature(Vector* P, int n, int dim, int type_descr, bool save_position=false)
{
    int i, j, init, end, init2, end2, dim_add;  
    double **features, *window, *coefs;

    double* p = P->get_elems();
    int N = P->get_n(); //vector length
    
    int w = round(N /(n*1.0)); // window
    double fact_add = 0.5;//ratio of window  for adding to both sides    
    // additional overlap
    int w_add = round(w * fact_add);
     
    // dim must not be greater than the size of window
    dim = std::min(dim, w ); 
    
    features=new double*[n];  
    
    init = 0;
    for(i=0;i<n;i++)
    {
        end=round(N*(i+1.0)/n);
        
        init2 = std::max(init - w_add, 0);
        end2  = std::min(end + w_add, N);

        window=new double[end2 - init2];
        for (j = init2; j < end2  ; j++)        
            window[j - init2] = p[j];
        
        features[i] = new double[dim];  
        
        dim_add =  save_position? 1 : 0;     
               
        coefs = apply_descriptor(window, end2 - init2, dim - dim_add, type_descr);
                        
        for(j = dim_add; j < dim ; j++)
        {
            features[i][j] = coefs[j - dim_add];                
        }

        init = end;

        if(save_position)
        {
            features[i][0] = equivalence(i, 0, n - 1, -1, 1);
        }       
                
        delete[] window;
        delete[] coefs;
    }
        
    return new Signature(features, NULL, n, dim, P->get_id_file());    
}


/**
 * extract n feature vectors from P
 */
std::vector<Object*> get_local_features(Vector* P, int n, int &dim,int type_descr, bool save_position=false)
{
    double *feature, *window, *coefs;
    int i, j, init, end, init2, end2, dim_add;  

    double* p = P->get_elems();
    int N = P->get_n(); //vector length

    int w = round(N /(n*1.0)); 
    double fact_add = 0.5;
    
    int w_add = round(w * fact_add);
        
    dim = std::min(dim, w ); 
    
    std::vector<Object*> features(n);
        
    init = 0;
    for(i = 0; i < n; i++)
    {
        end=round(N*(i+1.0)/n);
        
        init2 = std::max(init - w_add, 0);
        end2  = std::min(end + w_add, N);

        window=new double[end2 - init2];
        for (j = init2; j < end2  ; j++)        
            window[j - init2] = p[j];
        
        feature = new double[dim];
        
        dim_add = save_position? 1 : 0;
               
        coefs   = apply_descriptor(window, end2 - init2, dim - dim_add, type_descr);
                        
        for(j = dim_add; j < dim ; j++)
        {
            feature[j] = coefs[j - dim_add];                
        }

        if(save_position)
        {
            feature[0] = equivalence(i, 0, n - 1, -1, 1);           
        }   

        features[i] = new Vector(feature, dim, i);  

        init = end;
                
        delete[] window;
        delete[] coefs;
    }

    return features;
}


////////////////////// Heuristics for getting Signatures from UTS ///////////////////////////

/* 
 * @param  P             [data vector]
 * @param  nw            [size of the sliding window]
 * @param  overlap       [size of the overlap]
 * @param  n             [number total of feature vectors] 
 * @param  dim           [dimension of feature space]
 * @param  type_descr    [type_descr of the descriptor]
 * @param  n_clusters    [number of clusters]
 * @param  save_position [if we want to save the referential position]
 * */ 


/** heuristic 1: get the centroids by kmeans-technique 
  * in order to keep the referencial position of each local feature (SQFD)
  */
Object* get_signature_kmeans(Vector* P, int n, int dim, int n_clusters, int type_descr, bool save_position=false)
{   
    Lineal *dist = new Lineal(2);

    // get local features 
    std::vector<Object*> features = get_local_features(P, n, dim, type_descr, save_position);

    /////////////////// clustering kmeans ///////////////////////       
    vector<Object*> centroids_init(n_clusters);
    for(int i = 0; i < n_clusters; ++i) {       
        // seleccionar los feature centrales como centroides
        centroids_init[i] = features[ (i + 0.5) * n / n_clusters   -  1];
    }   
    vector< vector<Object*> > clusters = KMeans::get_kmeans(features, centroids_init, dist, 15);

    ///// generate signature
    double **centroids = new double*[n_clusters];
    double* weights = new double[n_clusters];
    for(int i = 0; i < n_clusters; ++i) {
        weights[i]  = clusters[i].size();
        centroids[i] = new double[dim];
        memcpy(centroids[i],  ((Vector*)centroids_init[i])->get_elems(), sizeof(double) * dim);
    }

    std::vector<Object*>().swap(features);

    return new Signature(centroids, weights, n_clusters, dim, P->get_id_file());   
}


void order_centroids(vector< vector<Object*> > &clusters, vector<Object*> &centroids)
{
    int n = centroids.size();

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if(centroids[i]->get_id_file() > centroids[j]->get_id_file())
            {               
                std::swap(centroids[i], centroids[j]);
                std::swap(clusters[i], clusters[j]);
            }
        }
    }
}

/** heuristic 2: get the centroids by kmediods-technique 
  * in order to keep the referencial position of each local feature (SQFD)
  * then ordering the result centroids
  */
Signature* get_signature_kmeans_ordered(Vector* P,  int n, int dim, int n_clusters,  int type_descr, bool save_position=false)
{
    Distance *dist = new Lineal(2);

    // get local features 
    std::vector<Object*> features = get_local_features(P, n, dim, type_descr, save_position);

    /////////////////// clustering kmeans ///////////////////////       
    vector<Object*> centroids_init(n_clusters);
    for(int i = 0; i < n_clusters; ++i) {       
        // seleccionar los feature centrales como centroides
        centroids_init[i] = features[ (i + 0.5) * n / n_clusters   -  1];
    }

    vector< vector<Object*> > clusters = KMedoids::get_kmedoids(features, centroids_init, dist);

    ///// order the centroids and the clusters  
    order_centroids(clusters, centroids_init);
    
    ///// generate signature
    double **centroids = new double*[n_clusters];
    double *weights = new double[n_clusters];
    for(int i = 0; i < n_clusters; ++i) {
        weights[i]  = clusters[i].size();
        centroids[i] = new double[dim];
        memcpy(centroids[i],  ((Vector*)centroids_init[i])->get_elems(), sizeof(double) * dim);
    }
        
    std::vector<Object*>().swap(features);

    return new Signature(centroids, weights, n_clusters, dim, P->get_id_file());   
}


/* heuristica 4- agrupar en dos niveles: se extrae n_cluster feature vector como centroides, luego n feature vectors se asignan mediante la funcion de distancia base */
Signature* get_signature_two_levels(Vector* P, int n, int dim, int n_clusters, int type_descr, bool save_position=false)
{
    Distance *dist = new Lineal(2);

    // get local features 
    std::vector<Object*> features = get_local_features(P, n, dim, type_descr, save_position);

    // get centroids features 
    std::vector<Object*> centroids_init = get_local_features(P, n_clusters, dim, type_descr, save_position);

    int i,j, pos_menor;
    double d,d_max;

    double* weights = new double[n_clusters];

    for (j = 0; j < n_clusters; ++j)
        weights[j] = 1;

    for(i = 0; i < n ; ++i) 
    {
        pos_menor = 0;
        d_max = dist->d(features[i], centroids_init[0]);        

        for (j = 1; j < n_clusters; ++j)
        {
            d = dist->d(features[i], centroids_init[j]);
            if(d < d_max)
            {
                pos_menor = j;
                d_max = d;
            }
        }

        weights[pos_menor] += 1;
    }

    double **centroids = new double*[n_clusters];
    for(i = 0; i < n_clusters; ++i) 
    {   
        centroids[i] = new double[dim];
        memcpy(centroids[i],  ((Vector*)centroids_init[i])->get_elems(), sizeof(double) * dim);
    }

    std::vector<Object*>().swap(features);

    return new Signature(centroids, weights, n_clusters, dim, P->get_id_file());
}


/* heuristica 3- agrupar por agregation: se fijan los centroides y el resto se agregan al mediode de la izquierda o derecha mediante la funcion de distancia base */
 // centroides fijos 
Signature* get_signature_aggregation_fixed(Vector* P, int n, int dim, int n_clusters, int type_descr, bool save_position=false)
{
    Distance *dist = new Lineal(2);

    // get local features 
    std::vector<Object*> features = get_local_features(P, n, dim, type_descr, save_position);
    
    double* weights = new double[n_clusters];
    vector<Object*> centroids_init(n_clusters);
    centroids_init[0] = features[0]; // el primer feature es el primer centroid
    weights[0] = 1;
    for(int i = 1; i < n_clusters; ++i) 
    {       
        // seleccionar los feature centrales como centroides
        centroids_init[i] = features[ 1.0 * i * n / (n_clusters - 1)   -  1];       
        weights[i] = 1;
    }

    ////////////// asignar pesos segun la cercania  de similitud  ///////////////////           
    int init, end; 
    Object *P1, *P2;
    for(int i = 0; i < n_clusters - 1; ++i) {
        P1 = centroids_init[i];
        P2 = centroids_init[i+1];
        init = P1->get_id_file() + 1;
        end  = P2->get_id_file() -  1;

        for(int j = init; j <= end; ++j){
            if(dist->d(P1, features[j])  < dist->d(P2, features[j]))
                weights[i] += 1;
            else
                weights[i+1] += 1;
        }   
    }
        
    double **centroids = new double*[n_clusters];
    for(int i = 0; i < n_clusters; ++i) {           
        centroids[i] = new double[dim];
        memcpy(centroids[i],  ((Vector*)centroids_init[i])->get_elems(), sizeof(double) * dim);
    }

    std::vector<Object*>().swap(features);

    return new Signature(centroids, weights, n_clusters, dim, P->get_id_file());   
}

// se obtiene el mean como centroide
Signature* get_signature_aggregation_mean(Vector* P, int n, int dim, int n_clusters, int type_descr, bool save_position=false)
{
    Lineal *dist = new Lineal(2);
    KMeans* kmeans = new KMeans(dist);

    // get local features 
    std::vector<Object*> features = get_local_features(P, n, dim, type_descr, save_position);
    
    double* weights = new double[n_clusters];
    vector<Object*> centroids_init(n_clusters);
    vector<Object*> cluster;
        
    ///////// seleccionar el feature mas representativo de un intervalo ////////    

    int i, j, init = 0, end = -1;
    for(i = 0; i < n_clusters ; ++i) 
    {       
        end =  n * (1.0 + i) / n_clusters - 1;

        for(j=init;j<=end;j++)      
            cluster.push_back(features[j]); 
        
        centroids_init[i] = kmeans->mean(cluster);
        weights[i] = 1;

        cluster.clear();
        
        init = end + 1;
    }

    ////////////// asignar pesos segun la cercania  de similitud  ///////////////////       
    for (j = 0; j < n_clusters; ++j)
        weights[j] = 1;
    
    Object *P1, *P2;

    init = 0.5 * n / n_clusters   -  1;
    P1 = centroids_init[0];

    for(i = 1; i < n_clusters ; ++i) 
    {   
        end = (i + 0.5) * n / n_clusters   -  1;    
        P2 = centroids_init[i];

        for(j = init; j <= end; ++j)
        {
            if(dist->d(P1, features[j])  < dist->d(P2, features[j]))
                weights[i-1] += 1;
            else
                weights[i] += 1;
        }

        init = end + 1;
        P1 = P2;
    }

        
    double **centroids = new double*[n_clusters];
    for(i = 0; i < n_clusters; ++i) 
    {   
        centroids[i] = new double[dim];
        memcpy(centroids[i],  ((Vector*)centroids_init[i])->get_elems(), sizeof(double) * dim);
    }

    std::vector<Object*>().swap(features);

    delete kmeans;

    return new Signature(centroids, weights, n_clusters, dim, P->get_id_file());   
}


// se obtiene el mediode como centroide
Signature* get_signature_aggregation_mediod(Vector* P, int n, int dim, int n_clusters, int type_descr, bool save_position=false)
{
    Distance *dist = new Lineal(2);
    KMedoids* kmediods = new KMedoids(dist);

    // get local features 
    std::vector<Object*> features = get_local_features(P, n, dim, type_descr, save_position);
    
    double* weights = new double[n_clusters];
    vector<Object*> centroids_init(n_clusters);
    vector<Object*> cluster;
        
    ///////// seleccionar el feature mas representativo de un intervalo ////////    

    int i, j, init = 0, end = -1;
    for(i = 0; i < n_clusters ; ++i) 
    {       
        end =  n * (1.0 + i) / n_clusters - 1;

        for(j=init;j<=end;j++)      
            cluster.push_back(features[j]); 
        
        centroids_init[i] = kmediods->medoide(cluster);
        weights[i] = 1;

        cluster.clear();
        
        init = end + 1;
    }

    ////////////// asignar pesos segun la cercania  de similitud  ///////////////////       
    Object *P1, *P2;
    for(i = 0; i < n_clusters - 1; ++i) 
    {
        P1 = centroids_init[i];
        P2 = centroids_init[i+1];
        init = P1->get_id_file() + 1;
        end  = P2->get_id_file() -  1;

        for(j = init; j <= end; ++j){
            if(dist->d(P1, features[j])  < dist->d(P2, features[j]))
                weights[i] += 1;
            else
                weights[i+1] += 1;
        }   
    }

    if(centroids_init[0]->get_id_file() > 0) // asignar los primeros features 
        weights[0] +=  centroids_init[0]->get_id_file();

    if(centroids_init[n_clusters - 1]->get_id_file() < (int) features.size() - 1) // asignar los ultimos features 
        weights[n_clusters - 1] +=  features.size() - centroids_init[n_clusters - 1]->get_id_file();

        
    double **centroids = new double*[n_clusters];
    for(i = 0; i < n_clusters; ++i) 
    {   
        centroids[i] = new double[dim];
        memcpy(centroids[i],  ((Vector*)centroids_init[i])->get_elems(), sizeof(double) * dim);
    }

    std::vector<Object*>().swap(features);

    delete kmediods;

    return new Signature(centroids, weights, n_clusters, dim, P->get_id_file());   
}


// se usa local feature para obtener los centroides
Signature* get_signature_aggregation_feature(Vector* P, int n, int dim, int n_clusters, int type_descr, bool save_position=false)
{
    Distance *dist = new Lineal(2);

    // get local features 
    std::vector<Object*> features = get_local_features(P, n, dim, type_descr, save_position);

    // get centroids features 
    std::vector<Object*> centroids_init = get_local_features(P, n_clusters, dim, type_descr, save_position);


    int i,j, init, end; 
    Object *P1, *P2;

    double* weights = new double[n_clusters];

    for (j = 0; j < n_clusters; ++j)
        weights[j] = 1;

    init = 0.5 * n / n_clusters   -  1;
    P1 = centroids_init[0];

    for(i = 1; i < n_clusters ; ++i) 
    {   
        end = (i + 0.5) * n / n_clusters   -  1;    
        P2 = centroids_init[i];

        for(j = init; j <= end; ++j)
        {
            if(dist->d(P1, features[j])  < dist->d(P2, features[j]))
                weights[i-1] += 1;
            else
                weights[i] += 1;
        }

        init = end + 1;
        P1 = P2;
    }

    double **centroids = new double*[n_clusters];
    for(i = 0; i < n_clusters; ++i) 
    {   
        centroids[i] = new double[dim];
        memcpy(centroids[i],  ((Vector*)centroids_init[i])->get_elems(), sizeof(double) * dim);
    }

    std::vector<Object*>().swap(features);

    return new Signature(centroids, weights, n_clusters, dim, P->get_id_file());   
}


////////////////////// Get Local Features & Signature from MTS //////////////////////////
 /* 
 * @param  Pinv          [multivariate time series]
 * @param  D             [number of variables]
 * @param  nw            [size of the sliding window]
 * @param  overlap       [size of the overlap]
 * @param  n             [number total of feature vectors] 
 * @param  dim           [dimension of feature space]
 * @param  type_descr    [type_descr of the descriptor]
 * @param  save_position [if we want to save the referential position]
 * */  

/**
 * select a binary simbol for correlation value [-1, 1]
 */
inline int select_symbol_corr(double corr)
{
    int val=1;
    if(corr >= 0.5) val = 3;
    else if(corr <= -0.5) val = 0;
    
    return val;
}

/**
 * get the weight of the relationship between the window components
 */
inline double apply_weight(double** window, int D, int N, int type_descr_weight=0)
{   
    double w = 1;    
  
    if(type_descr_weight == 1) // correlation
    {       
        w=0; 
        for (int i = 0; i < D - 1; ++i)
        {
            w += correlation(window[i], window[i+1], N);
        }
        w = 1 + w/(D-1);    
    }  
    else if(type_descr_weight == 2) // correlation & binary (cost based on hamming distance)
    {
        unsigned long wb = 0;
        for (int i = 0; i < D - 1; ++i)
        {
            w = correlation(window[i], window[i+1], N);
            wb = wb << 2;
            wb += select_symbol_corr(w);            
        }
        w = wb;     
    }
    
    return w;   
}


/** Get the signature receiving as parameter the sliding windows and the overlap */
Signature* get_signature(Matrix* Pinv, int nw, int overlap, int dim,  int type_descr, int type_weight)
{
    double  **window;
    int i, j, k, init, end, new_dim = 0; 

    double** ma = Pinv->get_elems();
    int D = Pinv->get_n();
    int N = Pinv->get_D();
    
    // number of feature vectors
    int n=ceil((N-nw)/(1.0*overlap))+1;

    // dim must not be greater than the size of window
    dim = std::min(dim, nw); 
        
    double **features = new double*[n];
    double* weights = new double[n];
        
    for(i = 0; i < n; i++)
    {
        init = i*overlap;
        end = std::min(init + nw, N);

        window = new double*[D];
        for (j = 0; j < D  ; j++)       
        {
            window[j]=new double[end - init]; 
            for (k = init; k < end  ; k++)      
                window[j][k - init] = ma[j][k];
        }

        if(type_descr == 0)
            features[i] = apply_stats(window, D, end - init, new_dim);
        else        
            features[i] = apply_descriptor(window, D, end - init, dim, new_dim, type_descr);


        weights[i] = apply_weight(window, D, end - init, type_weight);
         
                
        delete[] window;   
    }

    return new Signature(features, weights, n, new_dim, Pinv->get_id_file());
}

/** Get the signature receiving as parameter the number of feature vectors */
Signature* get_signature(Matrix* Pinv, int n, int dim, int type_descr, int type_weight)
{
    double  **window;
    int i, j, k, init, end, init2, end2, new_dim = 0; 

    double** ma = Pinv->get_elems();
    int D = Pinv->get_n();
    int N = Pinv->get_D();
    
    int nw = round(N /(n*1.0)); // window
    double fact_add = 0.5;//ratio of window  for adding to both sides    
    // additional overlap
    int overlap = round(nw * fact_add);
        
    double **features = new double*[n];
    double* weights = new double[n];
        
    init = 0;
    for(i = 0; i < n; i++)
    {
        end=round(N*(i+1.0)/n);
        
        init2 = std::max(init - overlap, 0);
        end2  = std::min(end + overlap, N);

        window = new double*[D];
        for (j = 0; j < D  ; j++)       
        {
            window[j]=new double[end2 - init2]; 
            for (k = init2; k < end2  ; k++)      
                window[j][k - init2] = ma[j][k];
        }

        if(type_descr == 0)
            features[i] = apply_stats(window, D, end2 - init2, new_dim);
        else        
            features[i] = apply_descriptor(window, D, end2 - init2, dim, new_dim, type_descr);


        weights[i] = apply_weight(window, D, end2 - init2, type_weight);
         
        init = end;

        delete[] window;   
    }

    return new Signature(features, weights, n, new_dim, Pinv->get_id_file());
}


/////////////// Hierarchical Feature from UTS //////////////////
/**
 * [Apply feature extraction hierarchically]
 * @param  P            [data vector]
 * @param  dim          [feature space  dimension]
 * @param  type_descr   [type_descr of descriptor]
 * @param  L      [resolution levels > 0]
 */
Signature* get_hie_features(Vector* P, int dim, int type_descr, int L = -1, bool by_levels=false)
{
    double *window, w;
    int n, init, end, i,j,l, k = 0;

    double* p = P->get_elems();
    int N = P->get_n(); //vector length


    int Lmax = log2(N/(1.0 * dim)  +  1); // max level of p    
    
    if(L > 0)
        L = std::min(Lmax, L); // adjust the user level  
    else
        L = Lmax;

    int total = pow(2, L) - 1; // total number of feature vectors
    
    double** features = new double*[total];
        
    for (l = 1; l <= L; ++l)
    {
        n = pow(2, l - 1);
        w = 1.0*N/n;
       
        window = new double[(int)ceil(w)];

        init = 0;
        for (i = 0; i < n; ++i)
        {
            end=std::min((int) round((i+1) * w), N);
        
            for (j = init; j < end; ++j)                           
            {
                window[j-init] = p[j]; 
            }    
            
            features[k]  = apply_descriptor(window, end - init, dim, type_descr);     
            
            init = end; 
            k++;           
        }

        delete[] window;
    }

    return new Signature(features, NULL, total, dim, P->get_id_file(), by_levels);
}


/////////////// Hierarchical Feature from MTS //////////////////
/**
 * [Apply feature extraction hierarchically]
 * @param  Pinv         [data matrix]
 * @param  dim          [feature space  dimension]
 * @param  type_descr   [type_descr of descriptor]
 * @param  L            [resolution levels > 0]
 */
Object* get_hie_features(Matrix* Pinv, int dim, int type_descr, int L = -1, bool by_levels=false)
{
    
    double **window, w;
    int n, init, end, i,j,l,d, new_dim = -1,  k = 0;

    double** ma = Pinv->get_elems();
    int D = Pinv->get_n();
    int N = Pinv->get_D(); 

    int Lmax = log2(N/(1.0 * dim)); // max level of p    
    
    if(L > 0)
        L = std::min(Lmax, L); // adjust the user level  
    else
        L = Lmax;

    int total = pow(2, L+1) - 1; // total number of feature vectors 
    double** features = new double*[total];
        
    for (l = 0; l <= L; ++l)
    {
        n = pow(2, l);
        w = 1.0*N/n;
       
        window = new double*[D];
        for (d = 0; d < D; d++)
            window[d] = new double[(int)ceil(w)];
        

        init = 0;
        for (i = 0; i < n; ++i)
        {
            end=std::min((int) round((i+1) * w), N);
        
            for (j = init; j < end; ++j)
                for (d = 0; d < D  ; d++)                            
                {
                    window[d][j-init] = ma[d][j]; 
                }    
            
            features[k]  = apply_descriptor(window, D, end - init, dim, new_dim, type_descr);  

            init = end; 
            k++;           
        }

        for (d = 0; d < D; d++)
            delete[] window[d];
        delete[] window;
    }

    return new Signature(features, NULL, total, new_dim, Pinv->get_id_file(), by_levels);
}

