#ifndef KMEANS_H
#define KMEANS_H

#include "../distances/lineal.h"
#include "../objects/object.h"
#include "../objects/vector.h"

#include <vector>

using namespace std;

class KMeans
{
private:
    //vectores caracteristicos
	vector<Object*> objects;
	
    // distancia lineal si es means, cualquier si es mediode
    Lineal* dist; // with symmetric property

public:
    KMeans(Lineal* _dist);
	KMeans(vector<Object*> _objects, Lineal* _dist);
    virtual ~KMeans();
    
    vector< vector<Object*> > run(vector<Object*> &centroids, int max_iter=25, double umbral=0.001);

    vector<Object*> init_centroids_rand(int n_cluster);

    vector<Object*> init_centroids_sort(int n_cluster);    
	
    Object* mean(vector<Object*> );

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
    
    double sum_dists(vector<Object*> objects, Object* query);
    int nearest(vector<Object*> objects, Object* query);

    
public:
    static vector< vector<Object*> > get_kmeans(vector<Object*> objects, vector<Object*> &centroids, Lineal* dist, int max_iter=25, double umbral=0.001);
};



#endif // KMEANS_H
