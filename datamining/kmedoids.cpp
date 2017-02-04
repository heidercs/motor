#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sstream>
#include <algorithm>

#include "kmedoids.h"
#include "../util/functions.h"
#include "../indexing/pairs.h"

KMedoids::KMedoids(vector<Object*> _objects, Distance* _dist)
{
	this->objects = _objects;
	this->dist = _dist;
}

KMedoids::KMedoids(Distance* _dist)
{	
	this->dist = _dist;
}

KMedoids::~KMedoids()
{
    if(cache != NULL)
	  free_cache();
}


vector< vector<Object*> > KMedoids::get_kmedoids(vector<Object*> objects, vector<Object*> &centroids, Distance* dist, int max_iter)
{
	KMedoids kmedoids(objects, dist);

    kmedoids.build_cache();

    return kmedoids.run_pam(centroids, max_iter);
}

void KMedoids::build_cache()
{
    int N = objects.size();

    this->original_index=new int[N];
    this->cache = new double*[N];

    for(int i=0; i<N; i++)
    {
    	//backup the original index_file
        original_index[i] = objects[i]->get_id_file();
        // set the new index file
        objects[i]->set_id_file(i);

        //initialize the cache
        cache[i] = new double[N];
    }

    // generating the cache
    double d;    
    for(int i=0; i<N; i++)
    {
        cache[i][i] = 0;
        for(int j=i+1; j<N; j++)
        {
            d = dist->d(objects[i], objects[j]);
            cache[i][j] = d; 
            cache[j][i] = d; 
        }
    }

}

void KMedoids::free_cache()
{
    if(cache == NULL) return;

    for(int i=0; i<(int)objects.size(); i++)
    {
        // recovery the original index
        objects[i]->set_id_file(original_index[i]);

        // free cache
        delete[] cache[i];
    }

    delete[] cache;
    delete[] original_index;
}

vector< vector<Object*> > KMedoids::associate_centroid(vector<Object*> &centroids)
{
	int n_cluster = centroids.size();
	int pos, N = (int) objects.size();

	vector< vector<Object*> > clusters(n_cluster);

	// agrupar los objects respecto a los centroids                 
    for(int i=0; i<N; i++)
    {
        // buscar el centroide mas cercano
        pos = nearest(centroids, objects[i]);

        //agregar a grupo seleccionado            
        clusters[pos].push_back(objects[i]);
    }

    // mezclar los clusters chicos
    merge_shorters(clusters, centroids, n_cluster, 1);

    return clusters;
}

vector< vector<Object*> > KMedoids::run_pam(vector<Object*> &centroids, int max_iter)
{
    vector< vector<Object*> > clusters;
    vector<Object*> new_centroids;    
    bool band;

    if((int)objects.size() < (int)centroids.size())  return clusters;//Validar que exista mas objects que numero de clusters

    double new_cost, prev_cost = DBL_MAX;

    clusters  = associate_centroid(centroids);
    prev_cost = WSS(clusters, centroids);

    int iter=0;
    do
    {
    	iter++;
        //std::cout<<"> Iter:"<<iter<<std::endl;
        
        // calcular nuevos centroides         
        for(int j=0; j<(int)centroids.size(); j++)
        {
        	//std::cout<<"\t- Grupo "<<j<<":"<<clusters[j].size()<<"\n";
			centroids[j] = medoide(clusters[j]);		          
        }

        clusters = associate_centroid(centroids);

        // verificar si los nuevos centroides generan un costo diferente
        new_cost = WSS(clusters, centroids);
        //std::cout<<"\t.. new cost:"<<new_cost<<std::endl;
        band = fabs(new_cost - prev_cost) >= 0.01;
        
        prev_cost = new_cost;

    }while(band && iter <= max_iter);

    return clusters;
}


void KMedoids::merge_shorters(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int n_cluster, int limit)
{			
	vector<Object*> cluster_merge;
	int pos;
	int i,j=0;   

    while(j< (int)clusters.size())
    {
    	if((int)clusters[j].size() <= limit)	
    	{
    		cluster_merge = clusters[j];  

			centroids.erase(centroids.begin() + j);
			clusters.erase(clusters.begin() + j);	
				    
			pos = -1;
			for (i = 0; i < (int)cluster_merge.size(); i++)
			{
			    pos = nearest(centroids, cluster_merge[i]);
				clusters[pos].push_back(cluster_merge[i]);
			}
    	}	
    	else j++;
    }

    if( (int)clusters.size() < n_cluster)
    {    
        for (i = 0; i < (int)centroids.size(); i++)
            centroids[i] = medoide(clusters[i]);

        // dividir los clusters grandes para completar el numero de clusters
        while((int)clusters.size() < n_cluster)
        {
        	pos = max_elements(clusters);	
        	split_cluster(clusters, centroids, pos); // dividir        	
        }	
    }
}

/**
 * [KMedoids::max_elements retorna la posicion del cluster con mas elementos] 
 */
int KMedoids::max_elements(vector< vector<Object*> > clusters)
{
	int count_max = 0, pos = -1;
	for(int i = 0; i < (int)clusters.size(); ++i) 
	{
		if((int) clusters[i].size() > count_max )
		{
			count_max = (int) clusters[i].size();
			pos = i;
		}
	}
	return pos;
}



/**
 * Inicializa el conjunto de centroides tomando
 * aleatoriamente N_CLUSTERS elementos
 */
vector<Object*> KMedoids::init_centroids_rand(int n_cluster, int seed)
{
    vector<Object*> centroids = vector<Object*>(n_cluster);

    if((int)objects.size() < n_cluster) return centroids;//validar tamanio de objects

    int desp=objects.size()/n_cluster;//Tamanio del desplazamiento
    int posRamdom=0;

    srand (time(NULL) * seed);    
    
    for(int i=0; i < n_cluster;  i++)
    {
        posRamdom=(int) (rand()%desp) ;        
        centroids[i]=objects[i*desp + posRamdom];
    }

    return centroids;
}


bool KMedoids::belongs(Object *p, vector<Object*> objs)
{
	bool result = false;
	for(int i=0; i < (int)objs.size() && !result; i++)
	{
		if(objs[i]->get_id_file() == p->get_id_file())
			result = true;
	} 
	
	return result;
}


/**
 * Inicializa el conjunto de centroides que maximisan la distancia entre si
 */
vector<Object*> KMedoids::init_centroids_dist(int n_cluster)
{
    vector<Object*> centroids;
    double d, d_max;
    int pos = -1, N = (int)objects.size();
    
    centroids = select_two_centroids(objects);
    
    while((int)centroids.size() < n_cluster)
    {
		d_max = 0;
		pos = -1;
		for(int i=0; i < N; i++)
		{
			if(belongs(objects[i], centroids)) continue;
			
            d = sum_dists(centroids,  objects[i]);

			if(d >= d_max)
			{
				d_max = d;
				pos = i;	
			}
		}
		centroids.push_back(objects[pos]);
	}   
    
    return centroids;
}


/**
 * Inicializa el conjunto de centroides mediante ordenacion
 */
vector<Object*> KMedoids::init_centroids_sort(int n_cluster)
{
    int i, j, N = (int)objects.size();

    double* summ = new double[N];
    for(i = 0; i < N; ++i)
    {
    	summ[i] = sum_dists(objects, objects[i]);
    }
    
    std::vector<PairObjectDistance> ranking(N);
    for (i= 0; i < N; ++i)
    {
    	ranking[i] = PairObjectDistance(objects[i], 0);
	    for (j = 0; j < N; ++j)
	    {
	    	ranking[i].dist += apply_dist(objects[i], objects[j]) / summ[j];
	    }
	}

	std::sort(ranking.begin(), ranking.end());

	vector<Object*> centroids;
	for (i = 0; i < n_cluster; ++i)
	{
		centroids.push_back(ranking[i].object);
	}    
    
    return centroids;
}


vector<Object*> KMedoids::select_two_centroids(vector<Object*> objects)
{
	vector<Object*> centroids;
	int pos = -1, N = (int)objects.size();
	double d, d_max;
	
	srand ( time(NULL) ); 
	Object* temp = objects[ rand()%N ];
	
	d_max=0;
	for(int i=0; i<N; i++)
	{
		d = apply_dist(temp, objects[i]);
		if(d > d_max)
		{
			d_max = d;
			pos = i;
		}
	}
	// insert the first centroid
	centroids.push_back(objects[pos]);
	
	d_max = 0;
	for(int i=0; i<N; i++)
	{
		d = apply_dist(centroids[0], objects[i]);
		if(d > d_max)
		{
			d_max = d;
			pos = i;
		}
	}	
	// insert the second centroid
	centroids.push_back(objects[pos]);
	
	return centroids;
}

double KMedoids::sum_dists(vector<Object*> objects, Object* query)
{
    double d = 0;
    for(int k=0; k < (int)objects.size();k++)
    {
        d += apply_dist(objects[k], query);
    }
    return d;
}

int KMedoids::nearest(vector<Object*> objects, Object* query)
{
    int posMenor=-1;
    double d, distmenor=DBL_MAX;

    for(int j=0; j<(int)objects.size(); j++)
    {
        d = apply_dist(objects[j], query);
        if(d < distmenor)
        {
            distmenor = d;
            posMenor = j;
        }
    }

    return posMenor;
}


Object* KMedoids::medoide(vector<Object*> objects)
{
    int N = (int)objects.size();
    double sse,sse_min=DBL_MAX;
    int pos=-1;

    for(int i=0;i<N;i++)
    {
        //sse = sum_dists(objects, objects[i]);  
        sse = SSE(objects, objects[i]);        
        
        if(sse < sse_min)
        {
            pos=i;
            sse_min=sse;
        }
    }

    return objects[pos];
}


double KMedoids::SSE(vector<Object*> cluster, Object* centroid)
{
	double sse = 0;
	for (int i = 0; i < (int)cluster.size() ; i++)
	{        
		sse += pow(apply_dist(cluster[i], centroid), 2);         
	}
	return sse;
}

double KMedoids::WSS(vector< vector<Object*> > clusters, vector<Object*> centroids)
{
    double sse=0;    

    for(int i=0;i<(int) centroids.size();i++)
    {        
        sse += SSE(clusters[i], centroids[i]);
    }

    return sse;
}

double KMedoids::BSS(vector< vector<Object*> > clusters, vector<Object*> centroids)
{
	Object* C=NULL;
    double sse=0; 
    int n_cluster = (int) centroids.size();   

    C = medoide(centroids);

    for(int i=0;i<n_cluster;i++)
    {        	
		sse += clusters[i].size() * pow(apply_dist(centroids[i], C), 2);	
    }

    return sse;
}

double KMedoids::silhouette(vector< vector<Object*> > clusters, vector<Object*>  centroids, int ind_cluster)
{
    double db,da,d, dmin=DBL_MAX;
    int ib=0;
    vector<Object*> cluster;        
   
	//identificar el cluster mas cercano	   
	for (int i = 0; i < (int)centroids.size() ; i++)
	{        
		if(i == ind_cluster) continue;

		d = apply_dist(centroids[i], centroids[ind_cluster]);
        if(d < dmin)
        {
        	dmin = d;
        	ib = i;
        }
	}
    

	//promediar la distancia a todos los elemntos del cluster ib		
    db = sum_dists(clusters[ib], centroids[ind_cluster]);
    db = db/clusters[ib].size();

	//promedia la distancia a todos los ementos en el mismo cluster		
    da = sum_dists(clusters[ind_cluster], centroids[ind_cluster]);

    da=da/(clusters[ind_cluster].size() - 1.0);//no se considera al medoide

	return (db-da)/std::max(db,da);	
    //return 1.0 - da/db;
}


double KMedoids::silhouette(vector< vector<Object*> > clusters, vector<Object*>  centroids)
{
    double sil = 0;
    
    for(int i=0;i<(int)clusters.size();i++)
    {
      sil += silhouette(clusters, centroids, i);
    }

    return sil;
}


bool KMedoids::pos_processing(vector< vector<Object*> > &clusters, vector<Object*> &centroids)
{
	int op, pos, n_cluster = (int)centroids.size();	
	vector<double> sse(n_cluster);
    double sse_max;
	
	for(int i=0;i < (int) centroids.size(); i++)
	{		
		sse[i] = SSE(clusters[i], centroids[i]);
	}	
	
	while(true)
	{

		////////////////////////////////////////////////////////////////////
		sse_max = 0;
        cout<<"Cluster\tsse\tsize\tavg_sse\n";
		for (int i = 0; i < n_cluster; i++)
		{
            cout<<i<<"\t"<<sse[i]<<"\t"<<clusters[i].size()<<"\t"<<sse[i]/clusters[i].size()<<endl;
            sse_max += sse[i]/clusters[i].size();
		}
		cout<<"\tSSE Total:"<<sse_max/n_cluster<<endl;        
		cout<<"Do you split or merge? (1 | 2):";
		cin>>op;
		
		if(op == 1)	
			cout<<"What cluster do you want to split? :";
		else if(op == 2)	
			cout<<"What cluster do you want to merge? :";	
		else
			break;
				
		cin >> pos;
		////////////////////////////////////////////////////////////////////
				
		if(pos >= 0 && pos < n_cluster)
		{		
			//split cluster
			if(op==1)			
			{
				cout<<"Split cluster "<<pos<<endl;
			
				split_cluster(clusters, centroids, pos);
				
				//generate SSE for each cluster
				sse[pos] = SSE(clusters[pos], centroids[pos]);		
				
				sse.push_back(SSE(clusters[n_cluster-1], centroids[n_cluster-1]));	
				
				n_cluster++;
			}
			else if(op == 2)
			{
				cout<<"Merge cluster "<<pos<<endl;
				
				sse.erase(sse.begin() + pos);
				
				pos = merge_cluster(clusters, centroids, pos);
				
				sse[pos] = SSE(clusters[pos], centroids[pos]);				
								
				n_cluster--;
				
			}
			else
				break;
		}		
	}
		
	vector<double>().swap(sse);	
	
	return true;
}

void KMedoids::split_cluster(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int ind_split)
{	
	vector<Object*> cluster_split, new_cluster, two_centroids;
	double d1, d2;

	cluster_split = clusters[ind_split];
	two_centroids = select_two_centroids(cluster_split);
	
	int i=0;	
	while(i < (int)cluster_split.size())
	{
		d1 = apply_dist(two_centroids[0], cluster_split[i]);
		d2 = apply_dist(two_centroids[1], cluster_split[i]);
		
		if(d2 < d1)
		{
			new_cluster.push_back(cluster_split[i]);
			cluster_split.erase(cluster_split.begin() + i);
		}
		else
			i++;
	}	
	
	//add the second cluster to clusters
	clusters[ind_split] = cluster_split;
	clusters.push_back(new_cluster);
	
	//add the second centroid to centroids
	centroids[ind_split] = two_centroids[0];
	centroids.push_back(two_centroids[1]);
}

int KMedoids::merge_cluster(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int ind_merge)
{			
	vector<Object*> cluster_merge = clusters[ind_merge];  
	
	centroids.erase(centroids.begin() + ind_merge);
	clusters.erase(clusters.begin() + ind_merge);	

		    
    int pos_min = -1;
	
	for (int i = 0; i < (int)cluster_merge.size(); i++)
	{
        pos_min = nearest(centroids, cluster_merge[i]);
		clusters[pos_min].push_back(cluster_merge[i]);
	}

    for (int i = 0; i < (int)centroids.size(); i++)
    {
        centroids[i] = medoide(clusters[i]);
    }
	
	return pos_min;
}

double KMedoids::apply_dist(Object* p, Object* q)
{
	double d = -1;

	// si el cache ha sido creado, buscar su valor
    if(cache != NULL)       
		d = cache[p->get_id_file()][q->get_id_file()];
	else  
		d = dist->d(p, q);
			
	return d;
}

