#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sstream>
#include <algorithm>

#include "kmeans.h"
#include "../util/functions.h"
#include "../indexing/pairs.h"

KMeans::KMeans(vector<Object*> _objects, Lineal* _dist)
{
	this->objects = _objects;
	this->dist = _dist;
}

KMeans::KMeans(Lineal* _dist)
{	
	this->dist = _dist;
}

KMeans::~KMeans()
{

}

vector< vector<Object*> > KMeans::get_kmeans(vector<Object*> objects, vector<Object*> &centroids, Lineal* dist, int max_iter, double umbral)
{
	KMeans kmeans(objects, dist);

	return kmeans.run(centroids, max_iter, umbral);
}

vector< vector<Object*> > KMeans::run(vector<Object*> &centroids, int max_iter, double umbral)
{
    vector< vector<Object*> > clusters;
    vector<Object*> new_centroids;    
    int i,j,iter, pos;
    double d;
    bool band;

    int N = (int) objects.size();
    int n_cluster = (int) centroids.size();

    if(N < n_cluster)  return clusters;//Validar que exista mas objects que numero de clusters

    new_centroids = vector<Object*>(n_cluster);

    for(j=0; j<n_cluster && !band; j++)

    iter=0;
    do
    {
    	iter++;
        //std::cout<<"> Iter:"<<iter<<std::endl;
        
        clusters.clear();
        clusters = vector< vector<Object*> >(n_cluster);
              
        // agrupar los objects respecto a los centroids                 
        for(i=0; i<N; i++)
        {
            // buscar el centroide mas cercano
            pos = nearest(centroids, objects[i]);

            //agregar a grupo seleccionado            
            clusters[pos].push_back(objects[i]);
        }

        // mezclar los clusters chicos
        merge_shorters(clusters, centroids, n_cluster, 1);
        
        // calcular nuevos centroides 
        for(j=0; j<n_cluster; j++)
        {
        	//std::cout<<"- N Grupo "<<j<<":"<<clusters[j].size()<<"\n";
			new_centroids[j] = mean(clusters[j]);		          
        }
        
        // verificar nuevos centroides: si la distancia menor entre centroids es mayor que el umbral continuar agrupando
        band=false;
        for(j=0; j<n_cluster && !band; j++)
        {			            
            d = dist->d(new_centroids[j], centroids[j]);
            band = d > umbral;  //seguir agrupando							         
        }		

        if(band && iter<=max_iter)
        {
            centroids.swap(new_centroids);
			n_cluster = (int) centroids.size();           
        }
               
    
    }while(band && iter<=max_iter);

    return clusters;
}


void KMeans::merge_shorters(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int n_cluster, int limit)
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
	    {
	        centroids[i] = mean(clusters[i]);
	    }


	    // dividir los clusters grandes para completar el numero de clusters
	    while((int)clusters.size() < n_cluster)
	    {
	    	pos = max_elements(clusters);	
	    	split_cluster(clusters, centroids, pos); // dividir        	
	    }
	}	
}

/**
 * [KMeans::max_elements retorna la posicion del cluster con mas elementos] 
 */
int KMeans::max_elements(vector< vector<Object*> > clusters)
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
vector<Object*> KMeans::init_centroids_rand(int n_cluster)
{
    vector<Object*> centroids = vector<Object*>(n_cluster);

    if((int)objects.size() < n_cluster) return centroids;//validar tamanio de objects

    int desp=objects.size()/n_cluster;//Tamanio del desplazamiento
    int posRamdom=0;

    srand ( time(NULL) );    
    
    for(int i=0; i < n_cluster;  i++)
    {
        posRamdom=(int) (rand()%desp) ;        
        centroids[i]=objects[i*desp + posRamdom];
    }

    return centroids;
}


/**
 * Inicializa el conjunto de centroides mediante ordenacion
 */
vector<Object*> KMeans::init_centroids_sort(int n_cluster)
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
	    	ranking[i].dist += dist->d(objects[i], objects[j]) / summ[j];
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


vector<Object*> KMeans::select_two_centroids(vector<Object*> objects)
{
	vector<Object*> centroids;
	int pos = -1, N = (int)objects.size();
	double d, d_max;
	
	srand ( time(NULL) ); 
	Object* temp = objects[ rand()%N ];
	
	d_max=0;
	for(int i=0; i<N; i++)
	{
		d = dist->d(temp, objects[i]);
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
		d = dist->d(centroids[0], objects[i]);
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

double KMeans::sum_dists(vector<Object*> objects, Object* query)
{
    double d = 0;
    for(int k=0; k < (int)objects.size();k++)
    {
        d += dist->d(objects[k], query);
    }
    return d;
}

int KMeans::nearest(vector<Object*> objects, Object* query)
{
    int posMenor=-1;
    double d, distmenor=DBL_MAX;

    for(int j=0; j<(int)objects.size(); j++)
    {        
        d = dist->d(objects[j], query);

        if(d < distmenor)
        {
            distmenor = d;
            posMenor = j;
        }
    }

    return posMenor;
}



Object* KMeans::mean(vector<Object*> objects)
{	
    int D = objects[0]->get_n();
    int N = (int)objects.size();
    
    double* med=new double[D];

    for(int j=0; j<D; j++)  med[j] = 0;

    for(int i=0; i < N;i++)
    {        
        for(int j=0; j<D; j++)
        {
            med[j] += ((Vector*) objects[i])->get(j);
        }
    }

    for(int j=0; j<D; j++)
        med[j] = med[j]/(1.0*N);

    return new Vector(med, D, -1);
}



double KMeans::SSE(vector<Object*> cluster, Object* centroid)
{
	double sse = 0;
	for (int i = 0; i < (int)cluster.size() ; i++)
	{        
		sse += pow(dist->d(cluster[i], centroid), 2);         
	}
	return sse;
}

double KMeans::WSS(vector< vector<Object*> > clusters, vector<Object*> centroids)
{
    double sse=0;    

    for(int i=0;i<(int) centroids.size();i++)
    {        
        sse += SSE(clusters[i], centroids[i]);
    }

    return sse;
}

double KMeans::BSS(vector< vector<Object*> > clusters, vector<Object*> centroids)
{
	Object* C=NULL;
    double sse=0; 
    int n_cluster = (int) centroids.size();   
    
    C = mean(centroids);		
        
    for(int i=0;i<n_cluster;i++)
    {        	
		sse += clusters[i].size() * pow(dist->d(centroids[i], C), 2);	
    }

    return sse;
}

double KMeans::silhouette(vector< vector<Object*> > clusters, vector<Object*>  centroids, int ind_cluster)
{
    double db,da,d, dmin=DBL_MAX;
    int ib=0;
    vector<Object*> cluster;        
   
	//identificar el cluster mas cercano	   
	for (int i = 0; i < (int)centroids.size() ; i++)
	{        
		if(i == ind_cluster) continue;

		d = dist->d(centroids[i], centroids[ind_cluster]);
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
    
    da=da/clusters[ind_cluster].size();

	return (db-da)/std::max(db,da);	
    //return 1.0 - da/db;
}


double KMeans::silhouette(vector< vector<Object*> > clusters, vector<Object*>  centroids)
{
    double sil = 0;
    
    for(int i=0;i<(int)clusters.size();i++)
    {
      sil += silhouette(clusters, centroids, i);
    }

    return sil;
}


bool KMeans::pos_processing(vector< vector<Object*> > &clusters, vector<Object*> &centroids)
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

void KMeans::split_cluster(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int ind_split)
{	
	vector<Object*> cluster_split, new_cluster, two_centroids;
	double d1, d2;

	cluster_split = clusters[ind_split];
	two_centroids = select_two_centroids(cluster_split);
	
	int i=0;	
	while(i < (int)cluster_split.size())
	{
		d1 = dist->d(two_centroids[0], cluster_split[i]);
		d2 = dist->d(two_centroids[1], cluster_split[i]);
		
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

int KMeans::merge_cluster(vector< vector<Object*> > &clusters, vector<Object*> &centroids, int ind_merge)
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
        centroids[i] = mean(clusters[i]);
    }

	
	return pos_min;
}

