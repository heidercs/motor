#ifndef KMEDOIDS_H
#define KMEDOIDS_H

#include "../distances/distance.h"
#include "../objects/object.h"
#include "../objects/vector.h"

#include <vector>

using namespace std;

class KMedoids
{
private:
    //vectores caracteristicos
	vector<Object*> objects;
	
    // distancia lineal si es means, cualquier si es mediode
    Distance* dist; // with symmetric property

	double** cache;

    int* original_index;
	
public:
    KMedoids(Distance* _dist);
	KMedoids(vector<Object*> _objects, Distance* _dist);
    virtual ~KMedoids();

    void build_cache();
    void free_cache();
    
    vector< vector<Object*> > run_pam(vector<Object*> &centroids, int max_iter=15);

    vector<Object*> init_centroids_rand(int n_cluster, int seed=1);
    vector<Object*> init_centroids_dist(int n_cluster);    
    vector<Object*> init_centroids_sort(int n_cluster);    
	
    Object* medoide(vector<Object*> );

    bool pos_processing(vector< vector<Object*> > &clusters, vector<Object*> &centroids);
	      
    // grado de cohesion en un cluster 
    double SSE(vector<Object*> cluster, Object* centroid);	
    
    // grado de cohesion en todos los clusters
    double WSS(vector< vector<Object*> > clusters, vector<Object*> centroids);
    
    // grado de separacion de los clusters, a nivel de centroides
    double BSS(vector< vector<Object*> > clusters, vector<Object*> centroids);
    
    // grado de separacion y cohesion de un cluster respecto a otro cluster mas cercano
    double silhouette(vector< vector<Object*> > clusters, vector<Object*>  centroids, int ind_cluster);
    
    // grado de separacion y cohesion de los clusters
    double silhouette(vector< vector<Object*> > clusters, vector<Object*>  centroids);

private:
    void split_cluster(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int ind_split);
    int  merge_cluster(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int ind_merge);
    void merge_shorters(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int n_cluster, int limit=1);
    vector<Object*> select_two_centroids(vector<Object*> );
    int max_elements(vector< vector<Object*> > clusters);
    
    double apply_dist(Object* p, Object* q);
    double sum_dists(vector<Object*> objects, Object* query);
    int nearest(vector<Object*> objects, Object* query);
    vector< vector<Object*> > associate_centroid(vector<Object*> &centroids);

    bool belongs(Object *p, vector<Object*> objs);
    
public:
    static vector< vector<Object*> > get_kmedoids(vector<Object*> objects, vector<Object*> &centroids, Distance* dist, int max_iter=15);
        
};



#endif // KMEANS_H
