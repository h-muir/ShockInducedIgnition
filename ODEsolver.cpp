//semi-analytic ODE solver for chemistry terms
//non-dimensionalised
//standalone system

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <assert.h>
using namespace std;

#define gamma 1.4   //heat capacity ratio for standard air
#define Q 1.24      //heat release parameter (non-dimensionalised)
#define Ea 14.99    //activation energy (non-dimensionalised)

typedef vector<double> Vector;

template<typename T>
void printVector(vector<T> vec){
    int len = vec.size();
    for(int i; i<len; ++i){
        cout << vec[i] << endl;
    }
}//vector-array print function templated to print for any vector type


const double epsilon = 1/Ea;
const double A = 52.8e5;
const double t0 = gamma/(A*exp(-Ea)*Q*(gamma-1));

Vector T_sol;   //temperature solution vector
Vector c_sol;   //concentration solution vector
Vector t_vec;   //time vector


double cnew(double cn, double Tn, double dt, double Tr, double cr){
    double beta = Q*gamma*cr/Tr;
    double cnew = cn*exp((-epsilon*Tr*dt/beta)*exp((1/(epsilon*Tr))*(1-1/Tn)));
    return cnew;
}

double Tnew(double cn, double dt, double Tr, double cr){
    //dt = dt/tr;
    double beta = Q*gamma*cr/Tr;
    double Tnew = 1 + beta*(1-cn);
    return Tnew;
}

void writetofile(Vector t_vec, Vector c_sol, Vector T_sol, double nT){
    ofstream file;
    file.open("./data/ODEsol.dat");
    for(int i = 0; i < nT; ++i){
        file << t_vec[i] << ' ';
        file << c_sol[i] << ' ';
        file << T_sol[i] << '\n';
    }
    file.close();
}


void writetoGNUplot(){
    cout << "compile final time solution with: \"chmod u+x ODEsol.gp\" the run ./ODEsol.gp \n";
    cout << "generated file ODEsol.svg located in ./visualisation/" << endl;  
    ofstream file;
    file.open("ODEsol.gp");
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal svg \n";
    file << "set out './visualisation/ODEsol.svg' \n";
    file << "set key off \n";
    file << "set grid ls 0 \n";
    file << "set multiplot layout 1,2 \n";
    file << "set title \"concentration\" \n";
    file << "set xlabel \'time (non-dimensional)\' \n";
    file << "plot 'ODEsol.dat' using 1:2 with linespoints pt 6 ps 0.5 lc rgb \'blue\'  \n";
    file << "set title \"Temperature\" \n";
    file << "plot 'ODEsol.dat' using 1:3 with linespoints pt 6 ps 0.5 lc rgb \'magenta\'  \n";
    file << "exit";

}



int main(){

    double c_new=1, T_new=1, t=0;
    double dt = 0.01;
    double c_final = 1e-5;

    T_sol.push_back(T_new);         //intial T
    c_sol.push_back(c_new);         //initial c
    t_vec.push_back(0);             //zero time

    double Tr = 1.0;
    double cr = 1.0;
    double tr = pow(Tr,2)/(cr*gamma)*exp(1/epsilon*(1/Tr-1));

    double i = 1;
    while(t < 2){
        t += dt;  
        c_new = cnew(c_new, T_new, dt/tr, Tr, cr);
        T_new = Tnew(c_new, dt/tr, Tr, cr);
        t_vec.push_back(t);
        c_sol.push_back(c_new);
        T_sol.push_back(T_new);
        i += 1;
    }

    double nT = i;
    cout << "final t = " << t << endl;
    cout << "iterations = " << nT << endl;
    cout << "tr = " << tr << endl;
    cout << "t0 = " << t0 << endl;
    cout << "max T = " << T_sol[nT-1] << endl;

    writetofile(t_vec, c_sol, T_sol, nT);
    writetoGNUplot();


    return 0;
}
