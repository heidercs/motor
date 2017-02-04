#ifndef READFILE_H
#define READFILE_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdlib.h>

using namespace  std;


void calc_dimension_archivo(const char* name_file, int &fil, int &col)
{

    ifstream input(name_file,ios::in);

    if(!input.good()) {
        cout<<"El fichero "<<name_file<<" no pudo abrirse\n";
        return ;
    }

    bool primeraVez=true;

    char buffer[32194];
    char *tokens;

    fil=0,col=0;

    while(!input.eof()){
        input.getline(buffer,32194);
        if(strlen(buffer)==0) continue;

        fil++;

        if(primeraVez){
            tokens=strtok(buffer,"  \t");
            tokens=strtok(NULL,"  \t");
            while(tokens != NULL) {
              tokens = strtok(NULL," \t");
              col++;
            }
            primeraVez=false;
        }

    }

    input.close();
}

double** leer_archivo(const char* name_file,int &fil, int &col){

    ifstream input(name_file,ios::in);

    if(!input.good()) {
        cout<<"El fichero "<<name_file<<" no pudo abrirse\n";
        return NULL;
    }

    double** matriz=new double*[fil];

    char buffer[32194];
    char *tokens;

    int i=0,j=0;

    while(!input.eof() && i<fil){
        input.getline(buffer,32194);
        if(strlen(buffer)==0) continue;

        matriz[i]=new double[col];

        j=0;
        tokens=strtok(buffer," \t");
        tokens=strtok(NULL,"  \t");
        while(tokens != NULL && j<col) {
            matriz[i][j]=atof(tokens);
            tokens = strtok(NULL, " \t");
            j++;
        }

        i++;
    }

    input.close();

    return matriz;
}

#endif // READFILE_H
