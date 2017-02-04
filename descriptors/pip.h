#ifndef PIP_H
#define PIP_H

#include <iostream>
#include <vector>

#include "../util/list.h"

// Perceptual Important Points
class Pip
{
public:
    int index;
    double value;

    Pip(){}
    Pip(int in, double v)
    {
        index = in;
        value = v;
    }
           
    friend std::ostream& operator <<(std::ostream& ,  Pip );
    friend bool operator >(Pip, Pip );
    friend bool operator <(Pip, Pip );

    static std::vector<Pip> get_all_pip(double* p, int N);
    static std::vector<Pip> get_rec_pip(double* p, int N, int new_n);
    static double* get_seg_pip(double* p, int N, int new_n);
    static double* get_rec_pip_double(double* p, int N, int new_n);

private:
    static int found_pip(double* time_serie, int left, int right);
    static void get_rec_pip(double* p, List<Pip>* l, int izq, int der, int n);    
};




#endif // PIP_H
