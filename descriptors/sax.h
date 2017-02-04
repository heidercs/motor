#ifndef SAX_H
#define SAX_H


#include "../objects/vector.h"

class Sax
{

public:    
    static const int CARD_DEFAULT=8;

private:    
    //numero de regiones breakpoints
    int card;
    //numero de bits por cada simbolo
    int tam_bits;   
    

    //parameters of the gaussian distribution
    double mean;
    double sd;


    //array of breakpoints
    double* break_points;
    //tabla para la distance entre dos simbolos SAX
    double** lookup_table;
    
    //puntos de corte. "p" normalizado  N(mean,sd)     
    void compute_breakpoints();

public:
    Sax() : Sax(CARD_DEFAULT, 0, 1) { };
    Sax(int _card) : Sax(_card, 0, 1) { };
    Sax(int _card, double _mean, double _sd);
    ~Sax();

    int* run_sax(double*p, int n, int new_n);
    int* run_sax_from_paa(double* paa, int new_n);

    bool* get_sax_bits(double*p, int n, int new_n);    

    // bit symbol size
    int get_tambits();
   
    //numero de regiones
    int get_card();
    
    //construye la tabla de distances de los simbolos SAX
    void construir_lookup();
    // return the cost between two sax values
    double cost_sax(int s1, int s2);
    //funcion de llamada
    double get_quantil(int p); 
   
    //convertir a bits la representacion sax numerica
    bool* convert_bits(int* sax, int n);

public:
    //calcular la representacion sax numerica
    static int* get_sax(double*p, int n, int new_n, int card=CARD_DEFAULT);    
    // sax a double
    static double* get_sax_double(double*p, int n, int new_n, int card=CARD_DEFAULT);
    
	// valores cuantiles de la distribucion normal, usados para obtener los breakpoints bajo el area de la curva 
	// se uso el lenguaje R con la funcion qnorm
    static constexpr double quantil_2[1]{0}; // qnorm(1/2)
	static constexpr double quantil_3[2]{-0.4307273, 0.4307273}; // qnorm(1/3), qnorm(2/3)
	static constexpr double quantil_4[3]{-0.6744898, 0, 0.6744898}; // qnorm(1/4), qnorm(2/4), qnorm(3/4)
	static constexpr double quantil_5[4]{-0.8416212, -0.2533471, 0.2533471, 0.8416212};
    static constexpr double quantil_6[5]{-0.9674216, -0.4307273, 0, 0.4307273, 0.9674216};
    static constexpr double quantil_7[6]{-1.0675705, -0.5659488, -0.1800124,  0.1800124,  0.5659488,  1.0675705};
    static constexpr double quantil_8[7]{-1.1503494, -0.6744898, -0.3186394,  0,  0.3186394,  0.6744898,  1.1503494};
	static constexpr double quantil_9[8]{-1.2206403, -0.7647097, -0.4307273, -0.1397103,  0.1397103,  0.4307273,  0.7647097,  1.2206403};
	static constexpr double quantil_10[9]{-1.2815516, -0.8416212, -0.5244005, -0.2533471,  0,  0.2533471,  0.5244005,  0.8416212,  1.2815516};
	static constexpr double quantil_16[15]{-1.53412054435255,-1.15034938037601,-0.887146559018876,-0.674489750196082,-0.488776411114669,-0.318639363964375,-0.157310684610171,0,0.157310684610171,0.318639363964375,0.488776411114669,0.674489750196082,0.887146559018876,1.15034938037601,1.53412054435255};	   

    ////////////////////////////////////////
    
};

class SaxFile
{
public:
    //representacion sax
    int* repr_sax;
    //tamaño de la representacion sax
    int n_sax;
    //serie orginal
    Vector* serie;
public:
    SaxFile();
    SaxFile(Vector* _serie, int* _repr_sax, int _n_sax);
    SaxFile(int* _repr_sax,int _n_sax);
    SaxFile(Vector* _serie);
};


#endif // SAX_H
