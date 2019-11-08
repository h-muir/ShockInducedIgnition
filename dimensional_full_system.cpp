//Dimensional full-system: HLLC solver for the Euler equations + semi-analytic ODE solver for chemistry
//This code is only intended for running the Toro tests in shockTests.H with chemistry: off
//full system with chemistry on should be run in "full_system.cpp"

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

const double A = 52.8e5;
const double t0 = gamma/(A*exp(-Ea)*Q*(gamma-1));
const double epsilon = 1/Ea;

typedef vector<double> Vector;

template<typename T>
void printVector(vector<T> vec){
    int len = vec.size();
    for(int i; i<len; ++i){
        cout << vec[i] << endl;
    }
}//vector-array print function templated to print for any vector type


//three class definitions for: 
//   - Primitive variables (W) 
//   - Conservative variables (U)
//   - Conservative fluxes (F) 

class Prim{
    public:
        Prim(double, double, double, double);
        double rho;
        double u;
        double p;
        double c;
};
Prim::Prim(double _rho, double _u, double _p, double _c){
    rho = _rho; u = _u; p=_p; c=_c;
}

class ConsU{
    public:
        ConsU(double, double, double, double);
        double rho;
        double rhou;
        double E;
        double rhoc;
};
ConsU::ConsU(double _rho, double _rhou, double _E, double _rhoc){
    rho = _rho; rhou = _rhou; E = _E; rhoc = _rhoc;
}

class ConsF{
    public:
        ConsF(double, double, double, double);
        double mass;
        double mom;
        double en;
        double cf;
};
ConsF::ConsF(double _mass, double _mom, double _en, double _cf){
    mass = _mass; mom = _mom; en = _en; cf = _cf;
}


//Print functions for classes-----------------------------------------------:
void printW(Prim W){
    printf("---------------\n W: \n rho: %f \n u:   %f \n p:   %f \n e:   %f \n c:    %f \n --------------\n",\
            W.rho, W.u, W.p, W.p/(W.rho*(gamma-1)), W.c);
}

void printU(ConsU U){
    printf("---------------\n U: \n rho:  %f \n rhou: %f \n E:    %f \n rhoc:    %f \n --------------\n"\
            ,U.rho, U.rhou, U.E, U.rhoc);
}

void printF(ConsF F){
    printf("---------------\n F: \n mass:  %f \n mom: %f \n en:    %f \n cf:    %f \n --------------\n"\
            ,F.mass, F.mom, F.en, F.cf);
}
//--------------------------------------------------------------------------/


//Conversion functions:

ConsU primtoconsU(Prim W){
    double rho = W.rho;
    double rhou = W.rho*W.u;
    double E = 0.5*W.rho*pow(W.u,2) + W.p/(gamma-1);
    double rhoc = W.rho*W.c;
    ConsU U(rho, rhou, E, rhoc);
    return U;
}

ConsF primtoconsF(Prim W){
    double mass = W.rho*W.u;
    double mom = W.rho*pow(W.u,2) + W.p;
    double E = 0.5*W.rho*pow(W.u,2) + W.p/(gamma-1);
    double en = W.u*(E + W.p);
    double cf = W.rho*W.u*W.c;
    ConsF F(mass, mom, en, cf);
    return F;
}

Prim consUtoprim(ConsU U){
    double rho = U.rho;
    double u = U.rhou/rho;
    double p = (gamma-1)*(U.E - 0.5*rho*pow(u,2));
    double c = U.rhoc/U.rho;
    Prim W(rho, u, p, c);
    return W;
}


#include "shockTests.H"

//for chemistry component:
Vector T_sol;   //temperature solution vector
Vector c_sol;   //concentration solution vector
Vector t_vec;   //time vector

double cnew(double cn, double Tn, double dt, double cr, double Tr){
    //double tr = pow(Tr,2)/(cr*gamma)*exp(1/epsilon*(1/Tr-1));
    //dt = dt/tr;  //scaling by tr places 'explosion' at t=1;
    double beta = Q*gamma*cr/Tr;
    double cnew = cn*exp((-epsilon*Tr*dt/beta)*exp((1/(epsilon*Tr))*(1-1/Tn)));
    return cnew;
}

double Tnew(double cn, double dt, double cr, double Tr){
    //double tr = pow(Tr,2)/(cr*gamma)*exp(1/epsilon*(1/Tr-1));
    //dt = dt/tr;  //scaling by tr places 'explosion' at t=1;
    double beta = Q*gamma*cr/Tr;
    double Tnew = 1 + beta*(1-cn);
    return Tnew;
}


ConsU delta(double omega, ConsU U_l, ConsU U_m, ConsU U_n, string scheme){
    /* omega = [-1, 1]
     * U_l = U_{i-1}, U_m = Ui, U_n = U{i+1}
     */
    ConsU delta_i_l(U_m.rho - U_l.rho, U_m.rhou - U_l.rhou, U_m.E - U_l.E, 
            U_m.rhoc - U_l.rhoc);
    ConsU delta_i_n(U_n.rho - U_m.rho, U_n.rhou - U_m.rhou, U_n.E - U_m.E,
            U_n.rhoc - U_m.rhoc);
    
    //alternative conventions for dividing by zero RHS slope:
    // -3 encoded to represent -ve infinity
    // +3 encoded to represent +ve infinity 
    double r_rho = delta_i_l.rho;
    if(delta_i_n.rho == 0){
        if(delta_i_l.rho == 0){
            r_rho = 0;           //both left and right slopes are zero
        }else if(delta_i_l.rho < 0){
            r_rho = -3;          //negative slope on left, zero slope on right
        }else{
            r_rho = 3;           //positive slope on left, zero slope on right
        }
    }else{r_rho = r_rho/delta_i_n.rho;}     //permissible division by delta_i_n
    double r_rhou =  delta_i_l.rhou;
    if(delta_i_n.rhou == 0){
       if(delta_i_l.rhou == 0){
            r_rhou = 0;          //both left and right slopes are zero
        }else if(delta_i_l.rhou < 0){
            r_rhou = -3;         //negative slope on left, zero slope on right
        }else{
            r_rhou = 3;          //positive slope on left, zero slope on right
        }
    }else{r_rhou = r_rhou/delta_i_n.rhou;}  //permissible division by delta_i_n
    double r_E = delta_i_l.E;
    if(delta_i_n.E == 0){ 
        if(delta_i_l.E == 0){
            r_E = 0;          //both left and right slopes are zero
        }else if(delta_i_l.E < 0){
            r_E = -3;         //negative slope on left, zero slope on right
        }else{
            r_E = 3;          //positive slope on left, zero slope on right
        }
    }else{r_E = r_E/delta_i_n.E;}           //permissible division by delta_i_n
    double r_rhoc = delta_i_l.rhoc;
    if(delta_i_n.rhoc ==0){
        if(delta_i_l.rhoc == 0){
            r_rhoc = 0;          //both left and right slopes are zero
        }else if(delta_i_l.rhoc < 0){
            r_rhoc = -3;         //negative slope on left, zero slope on right
        }else{
            r_rhoc = 3;          //positive slope on left, zero slope on right
        }
    }else{r_rhoc = r_rhoc/delta_i_n.rhoc;} 

    ConsU delta_i(0.5*(1+omega)*delta_i_l.rho+0.5*(1-omega)*delta_i_n.rho,
                  0.5*(1+omega)*delta_i_l.rhou+0.5*(1-omega)*delta_i_n.rhou,
                  0.5*(1+omega)*delta_i_l.E+0.5*(1-omega)*delta_i_n.E,
                  0.5*(1+omega)*delta_i_l.rhoc+0.5*(1-omega)*delta_i_n.rhoc);

    ConsU sigma_R(2/(1-omega+(1+omega)*r_rho), 2/(1-omega+(1+omega)*r_rhou),
                  2/(1-omega+(1+omega)*r_E), 2/(1-omega+(1+omega)*r_rhoc));

        //adjustments for infinity encoded cases:
    if(r_rho == -3 || r_rho == 3){
        sigma_R.rho = 0;
    }
    if(r_rhou == -3 || r_rhou == 3){
        sigma_R.rhou = 0;
    }
    if(r_E == -3 || r_E == 3){
        sigma_R.E = 0;
    }
    if(r_rhoc == -3 || r_rhoc == 3){
        sigma_R.rhoc = 0;
    } 

    if(scheme == "SuperBee"){
        ConsU delta_bar(0,0,0,0);
        //cout << "r_rho: " << r_rho << endl;
        if(r_rho < 0){
            delta_bar.rho = 0;
        }else if(0 <= r_rho && r_rho < 0.5){
            delta_bar.rho = 2*r_rho*delta_i.rho;
        }else if(0.5 <= r_rho && r_rho < 1){
            delta_bar.rho = 1*delta_i.rho;
        }else if(1 <= r_rho && r_rho < 2){
            delta_bar.rho = ((r_rho < sigma_R.rho)? 
                    r_rho*delta_i.rho : sigma_R.rho*delta_i.rho);
        }else{
            delta_bar.rho = ((sigma_R.rho < 2)? 
                    sigma_R.rho*delta_i.rho : 2*delta_i.rho);
;
        }
        if(r_rhou < 0){
            delta_bar.rhou = 0;
        }else if(0 <= r_rhou && r_rhou < 0.5){
            delta_bar.rhou = 2*r_rhou*delta_i.rhou;
        }else if(0.5 <= r_rhou && r_rhou < 1){
            delta_bar.rhou = 1*delta_i.rhou;
        }else if(1 <= r_rhou && r_rhou < 2){
            delta_bar.rhou = ((r_rhou < sigma_R.rhou)? 
                    r_rhou*delta_i.rhou : sigma_R.rhou*delta_i.rhou);
        }else{
            delta_bar.rhou = ((sigma_R.rhou < 2)? 
                    sigma_R.rhou*delta_i.rhou : 2*delta_i.rhou);
        }
        if(r_E < 0){
            delta_bar.E = 0;
        }else if(0 <= r_E && r_E < 0.5){
            delta_bar.E = 2*r_E*delta_i.E;
        }else if(0.5 <= r_E && r_E < 1){
            delta_bar.E = 1*delta_i.E;
        }else if(1 <= r_E && r_E < 2){
            delta_bar.E = ((r_E < sigma_R.E)? 
                    r_E*delta_i.E : sigma_R.E*delta_i.E);
        }else{
            delta_bar.E = ((sigma_R.E < 2)? 
                    sigma_R.E*delta_i.E : 2*delta_i.E);
        }

        if(r_rhoc < 0){
            delta_bar.rhoc = 0;
        }else if(0 <= r_rhoc && r_rhoc < 0.5){
            delta_bar.rhoc = 2*r_rhoc*delta_i.rhoc;
        }else if(0.5 <= r_rhoc && r_rhoc < 1){
            delta_bar.rhoc = 1*delta_i.rhoc;
        }else if(1 <= r_rhoc && r_rhoc < 2){
            delta_bar.rhoc = ((r_rhoc < sigma_R.rhoc)? 
                    r_rhoc*delta_i.rhoc : sigma_R.rhoc*delta_i.rhoc);
        }else{
            delta_bar.rhoc = ((sigma_R.rhoc < 2)? 
                    sigma_R.rhoc*delta_i.rhoc : 2*delta_i.rhoc);
        }


        return delta_bar;

    }else{

        return delta_i;

    }
}


class Cell{
    public:
        double aL, aR;                  //left and right sound speeds
        double SL, SR, Sstar, Splus;    //wave speeds
        ConsU UL, UR, ULstar, URstar;   //4 cons. states
        ConsF FL, FR, FLstar, FRstar;   //4 cons. fluxes
        ConsU U;                        //tbd state
        ConsF Fout;                     //flux out (F_i+1/2)
        
        //must declare constructor here:
        Cell(Prim WL, Prim WR):
            //must initialize all class objects here:
            UL(0,0,0,0), UR(0,0,0,0), ULstar(0,0,0,0), URstar(0,0,0,0),
            FL(0,0,0,0), FR(0,0,0,0), FLstar(0,0,0,0), FRstar(0,0,0,0),
            U(0,0,0,0), Fout(0,0,0,0) {
            //direct wave speed estimates:
            aL = pow(fabs(gamma*WL.p/WL.rho), 0.5); //non-physical states are permitted... 
            aR = pow(fabs(gamma*WR.p/WR.rho), 0.5); //...as an intermediate step...
            //... therefore fabs corrects for negative pressure and/or density
            SL = ((WL.u - aL < WR.u - aR)? WL.u - aL : WR.u - aR);
            SR = ((WL.u + aL > WR.u + aR)? WL.u + aL : WR.u + aR);
            Sstar = (WR.p - WL.p + WL.rho*WL.u*(SL-WL.u) - WR.rho*WR.u*(SR-WR.u))/\
                    (WL.rho*(SL-WL.u)-WR.rho*(SR-WR.u));
            Splus = ((fabs(WL.u)+aL > fabs(WR.u)+aR)? fabs(WL.u)+aL : fabs(WR.u)+aR);
            //W-->U and W-->F conversions:
            UL = primtoconsU(WL);
            UR = primtoconsU(WR);
            FL = primtoconsF(WL);
            FR = primtoconsF(WR);
            //Star region calculations:
            ULstar.rho = WL.rho*((SL-WL.u)/(SL-Sstar));
            ULstar.rhou = ULstar.rho*Sstar;
            ULstar.E = ULstar.rho*(UL.E/WL.rho+(Sstar-WL.u)*(Sstar + WL.p/(WL.rho*(SL-WL.u))));
            ULstar.rhoc = ULstar.rho*WL.c;
            //assert(ULstar.E >=0); //check that total energy is non-negative
            URstar.rho = WR.rho*((SR-WR.u)/(SR-Sstar));
            URstar.rhou = URstar.rho*Sstar;
            URstar.E = URstar.rho*(UR.E/WR.rho+(Sstar-WR.u)*(Sstar + WR.p/(WR.rho*(SR-WR.u))));
            URstar.rhoc = URstar.rho*WR.c;
            //assert(URstar.E >=0); //check that total energy is non-negative
            FLstar.mass = FL.mass + SL*(ULstar.rho - UL.rho);
            FLstar.mom = FL.mom + SL*(ULstar.rhou - UL.rhou);
            FLstar.en = FL.en + SL*(ULstar.E - UL.E);
            FLstar.cf = FLstar.mass*WL.c;
            FRstar.mass = FR.mass + SR*(URstar.rho - UR.rho);
            FRstar.mom = FR.mom + SR*(URstar.rhou - UR.rhou);
            FRstar.en = FR.en + SR*(URstar.E - UR.E);
            FRstar.cf = FRstar.mass*WR.c;
            }

};

void updateBoundary(vector<Prim> &W_vec, int n, string condition){
    W_vec[1] = W_vec[2];
    W_vec[0] = W_vec[1];
    W_vec[n-2] = W_vec[n-3];
    W_vec[n-1] = W_vec[n-2];

    if(condition == "left reflective wall"){
        W_vec[0].u = -W_vec[3].u;   //reflected shock BCs
        W_vec[1].u = -W_vec[2].u;   //reflected shock BCs        
    }else if(condition == "transmissive"){
        //transmissive conditions defined above
    }else{ cout << "!!! non-existent BC (typo?) " << endl; }
    
}

vector<Prim> updateW(vector<Prim> W_vec, int n, double &dt, double dx, double CFL, 
        double &Smax, bool MUSCL, string scheme, string BC){

    //MUSCL-Hancock parameters:
    ConsU UL(0,0,0,0);            
    ConsU Ul(0,0,0,0);
    ConsU Ui(0,0,0,0);
    ConsU Un(0,0,0,0);
    ConsU UR(0,0,0,0);
    ConsU UbarL(0,0,0,0);
    ConsU UbarR(0,0,0,0);
    Prim WbarL(0,0,0,0);
    Prim WbarR(0,0,0,0);
    ConsF FbarL(0,0,0,0);
    ConsF FbarR(0,0,0,0);    
    double omega = 0;
    ConsU delta_i(0,0,0,0);
    vector<Prim> WbarL_vec = W_vec;
    vector<Prim> WbarR_vec = W_vec; 
    
    //update paramters 
    vector<Prim> W_vec_new = W_vec; //initialised variable for subsequent updating  
    Cell cell0(W_vec[0], W_vec[0]);
    vector<Cell> cell_vec(n, cell0);
    ConsU Unew(0,0,0,0);
    Prim Wnew(0,0,0,0);

    //updateBoundary(W_vec, n, "left reflective wall");
    for(int i = 2; i<(n-2); ++i){
        //MUSCL-Hancock steps 1+2:
        Ul = primtoconsU(W_vec[i-1]); //left-side boundary
        Ui = primtoconsU(W_vec[i]);
        Un = primtoconsU(W_vec[i+1]); //right-side boundary
        
        
        if(MUSCL == 1){        
            delta_i = delta(omega, Ul, Ui, Un, scheme);
        }
        

        UL = ConsU(Ui.rho - 0.5*delta_i.rho,
                   Ui.rhou - 0.5*delta_i.rhou,
                   Ui.E - 0.5*delta_i.E,
                   Ui.rhoc - 0.5*delta_i.rhoc);
        UR = ConsU(Ui.rho + 0.5*delta_i.rho,
                   Ui.rhou + 0.5*delta_i.rhou,
                   Ui.E + 0.5*delta_i.E,
                   Ui.rhoc + 0.5*delta_i.rhoc);
        WbarL = consUtoprim(UL);
        WbarR = consUtoprim(UR);
        FbarL = primtoconsF(WbarL);
        FbarR = primtoconsF(WbarR);
        UbarL.rho = UL.rho + 0.5*(dt/dx)*(FbarL.mass - FbarR.mass);
        UbarL.rhou = UL.rhou + 0.5*(dt/dx)*(FbarL.mom - FbarR.mom);
        UbarL.E = UL.E + 0.5*(dt/dx)*(FbarL.en - FbarR.en);
        UbarL.rhoc = UL.rhoc + 0.5*(dt/dx)*(FbarL.cf - FbarR.cf);
        UbarR.rho = UR.rho + 0.5*(dt/dx)*(FbarL.mass - FbarR.mass);
        UbarR.rhou = UR.rhou + 0.5*(dt/dx)*(FbarL.mom - FbarR.mom);
        UbarR.E = UR.E + 0.5*(dt/dx)*(FbarL.en - FbarR.en);
        UbarR.rhoc = UR.rhoc + 0.5*(dt/dx)*(FbarL.cf - FbarR.cf);

        WbarL = consUtoprim(UbarL);
        WbarR = consUtoprim(UbarR);
        
        WbarL_vec[i] = WbarL;
        WbarR_vec[i] = WbarR;
        

    }
    for(int i = 0; i<n; ++i){
        //cout << i << endl; 
        if(i==n-1){
            cell_vec[i] = Cell(WbarR_vec[i], WbarR_vec[i]);
        }else{
            cell_vec[i] = Cell(WbarR_vec[i], WbarL_vec[i+1]);                
        }

        if(cell_vec[i].Splus > Smax){
            Smax = cell_vec[i].Splus;
        }
    }
    assert(Smax > 0);
    dt = CFL*dx/Smax;
    assert(dt > 0);     //check this is true rather than enforcing fabs(dt)
    Cell cell = cell0;  //any initialisation
    Cell prev_cell = cell0;  

    for(int i = 0; i<n; ++i){

        cell = cell_vec[i];
        Ui = primtoconsU(W_vec[i]);

        if( 0 < cell.SL){
            cell.U = cell.UL;
            cell.Fout = cell.FL;
        }else if(cell.SL <= 0 && 0 < cell.Sstar){
            cell.U = cell.ULstar;
            cell.Fout = cell.FLstar;
        }else if(cell.Sstar <= 0 && 0 < cell.SR){
            cell.U = cell.URstar;
            cell.Fout = cell.FRstar;
        }else if(cell.SR <= 0){
            cell.U = cell.UR;
            cell.Fout = cell.FR;
        }
        
        
        if(i == 0){
            prev_cell = cell;
        }
         
        //conservative update formula
        Unew.rho = Ui.rho + dt/dx * (prev_cell.Fout.mass - cell.Fout.mass);
        Unew.rhou = Ui.rhou + dt/dx * (prev_cell.Fout.mom - cell.Fout.mom);
        Unew.E = Ui.E + dt/dx * (prev_cell.Fout.en - cell.Fout.en);
        Unew.rhoc = Ui.rhoc + dt/dx * (prev_cell.Fout.cf - cell.Fout.cf);
        
        //variable conversion
        Wnew = consUtoprim(Unew);
        W_vec_new[i] = Wnew;
        prev_cell = cell; 
                    
    }

    updateBoundary(W_vec_new, n, BC);

    return W_vec_new;
}

vector< vector<Prim> > generate_solution(TestCase CASE, Vector &time_vec, 
        bool MUSCL, string scheme, string BC, bool wall_correction, bool chemistry, int ns){

    int n = CASE.n;
    double CFL = CASE.CFL;
    double T0 = CASE.T0;
    double Tf = CASE.Tf;
    double dx = CASE.dx;
    double T = T0;
    double L = CASE.L;
    double x0 = CASE.x0;
    Prim WL0 = CASE.WL;
    Prim WR0 = CASE.WR;
    Vector space = CASE.space;

    Cell cell0(WL0, WL0);
    vector<Cell> cell_vec(n, cell0);
    double Smax=0;
    double dt = 0;
    double dt_chem;                 //timestep for chemistry sub-cycling
    double t_chem;                  //time tracking for chemistry sub-cycle
    double c_new, T_new;            //chem
    Vector C_vec(n,1);
    Vector T_vec(n,1);
    //Vector C_vec_new = C_vec; 
    double p_sub;
    double cr, Tr, tr;

    vector< vector<Prim> > MATRIX;  //solution is generated as a space*time matrix 
                                    //of the primitive variables 
    vector<Prim> W_vec(n, WR0);     //spatial vector can be initialised with known n states
    vector<Prim> W_vec_new = W_vec; //initialised variable for subsequent updating 
    //adjusting left-side state (WL0) initialisation:
    for(int i=0; i < (int)(x0/L*n); ++i){ 
        //where x0 is the discontinuity point in the global Riemann problem
        W_vec[i] = WL0;
    }

    updateBoundary(W_vec, n, BC);

    MATRIX.push_back(W_vec);
    for(int i = 1; i < n; ++i){
        cell0 = Cell(W_vec[i-1], W_vec[i]);
        if(cell0.Splus > Smax){
            Smax = cell0.Splus;
        }
    }
    
    cout << "Smax0 = " << Smax << endl;
    assert(Smax > 0);
    dt = CFL*dx/Smax;

    time_vec.push_back(0);
    int nT = 1; //zero time is the first time step

    double tol = 1e-3;
    const double T0r = WL0.p/WL0.rho;
    bool shock_detection_flag = 0;
    bool wall_correction_flag = 0;
    bool correction_completed = 0;

    
    while(T < Tf){
        Smax = 0;
        W_vec_new = updateW(W_vec, n, dt, dx, CFL, Smax, MUSCL, scheme, BC);
        W_vec = W_vec_new;
        
        dt_chem = dt/ns;       //division by number of sub-cycle steps
        
        if(wall_correction){
            //wall heating correction:
            if(fabs((W_vec[10].rho - W_vec[8].rho)) > 0.2){
                //detection of sharp gradient ahead of critical cells
                shock_detection_flag = 1;
                //cout << "shock detected" << endl;
                //once the steep gradient is detected, the shock has been flagged to have passed
                //and the flag stays = 1 herein
            }
            if(shock_detection_flag && fabs((W_vec[10].rho - W_vec[8].rho)) < 0.002){
                //the cells upstream of the erroneous cells requiring correction
                //must have settled to the true value
                //when the two cells upstream of the error cells have a gradient < 0.2% 
                //the wall correction flag is switched on
                wall_correction_flag = 1;
                //cout << "wall correction flagged" << endl;
            }
            if(shock_detection_flag && wall_correction_flag && correction_completed == 0){
                //once the two conditions are triggered:
                //  (1) shock wave has passed
                //  (2) cells behind shock have settled
                //the wall correction is performed for density and concentration
                for(int i = 9; i >= 0; --i){
                    W_vec[i].rho = 1; 
                }
                cout << "correction completed" << endl;
                correction_completed = 1;
            }
        }

        //CHEMISTRY subcycle:
        if(chemistry){
            //chemistry: on/off
            for(int i = 0; i<n; ++i){
                //chemistry is only added to the system after the wall fix has occurred
                if(correction_completed){ 
                    t_chem = 0;
                    cr = W_vec[i].c; //C_vec[i];
                    Tr = W_vec[i].p/(W_vec[i].rho);
                    T_new = 1;
                    c_new = 1;
                    tr = pow(Tr,2)/(cr*gamma)*exp(1/epsilon*(1/Tr-1));
                    while(t_chem < dt){
                        t_chem += dt_chem;      //timestep is exactly matched with PDE dt
                        c_new = cnew(c_new, T_new, dt_chem/tr, cr, Tr);
                        T_new = Tnew(c_new, dt_chem/tr, cr, Tr);

                    }
                    T_vec[i] = T_new*Tr;
                    C_vec[i] = c_new*cr;
                    p_sub = W_vec[i].rho*T_vec[i];
                    W_vec[i].p = p_sub; 
                    W_vec[i].c = C_vec[i];                      
                }
            }
        }

        T += dt;
        cout << "T = " << T << endl;
        //current simulation time printed to terminal
        time_vec.push_back(T);        
        nT += 1; //time step counter        
        MATRIX.push_back(W_vec);            
    }
    
    cout << "(actual final time = " << T << "s)" << endl; 
    
    return MATRIX;
}

//(beginning) DATA AND VISUALISATION -----------------------------------------------------//

Vector writetofile(TestCase CASE, vector< vector<Prim> > MATRIX){
    Vector output(10); //returning: {u_min, u_max, rho_min, rho_max,
                                // p_min, p_max, e_min, e_max, c_min, c_max}
                                // for use in automatic plot scaling later
    Vector space = CASE.space;
    string name = CASE.datfile;
    int n = CASE.n;     //number of spatial nodes
    int nT = CASE.nT;   //number of time steps
    assert(MATRIX.size() == nT);
    assert(MATRIX[0].size() == n);
    ofstream ufile, rhofile, pfile, efile, cfile, finalT;
    ufile.open("./data/u_" + name);
    rhofile.open("./data/rho_" + name);
    pfile.open("./data/p_" + name);
    efile.open("./data/e_" + name);
    cfile.open("./data/c_" + name);
    finalT.open("./data/finalT" + name);
    double u_min=CASE.WL.u, u_max=u_min, rho_min=CASE.WL.rho, rho_max=rho_min; 
    double p_min=CASE.WL.p, p_max=p_min, e_min=p_min/(rho_min*(gamma-1)), e_max=e_min;
    double c_min=CASE.WL.c, c_max=c_min;
    //just initialising these to possible values
    double e;
    for(int i = 0; i < n; ++i){
        //first column of each file is the space vector
        ufile << space[i] << ' ';
        rhofile << space[i] << ' ';
        pfile << space[i] << ' ';
        efile << space[i] << ' ';
        cfile << space[i] << ' ';
        finalT << space[i] << ' ';
        for(int j = 0; j < nT; ++j){
            ufile << MATRIX[j][i].u << ' ';
            if(MATRIX[j][i].u < u_min){
                u_min = MATRIX[j][i].u;
            }
            if(MATRIX[j][i].u > u_max){
                u_max = MATRIX[j][i].u;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].u << ' ';
            }
            rhofile << MATRIX[j][i].rho << ' ';
            if(MATRIX[j][i].rho < rho_min){
                rho_min = MATRIX[j][i].rho;
            }
            if(MATRIX[j][i].rho > rho_max){
                rho_max = MATRIX[j][i].rho;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].rho << ' ';
            }

            pfile << MATRIX[j][i].p << ' ';
            if(MATRIX[j][i].p < p_min){
                p_min = MATRIX[j][i].p;
            }
            if(MATRIX[j][i].p > p_max){
                p_max = MATRIX[j][i].p;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].p << ' ';
            }

            e = MATRIX[j][i].p/(MATRIX[j][i].rho*(gamma-1));
            efile << e << ' ';
            if(e<e_min){
                e_min = e;
            }
            if(e>e_max){
                e_max = e;
            }
            if(j== nT-1){
                finalT << e << ' ';
            }
            cfile << MATRIX[j][i].c << ' ';
            if(MATRIX[j][i].c < c_min){
                c_min = MATRIX[j][i].c;
            }
            if(MATRIX[j][i].c > c_max){
                c_max = MATRIX[j][i].c;
            }
            if(j== nT-1){
                finalT << MATRIX[j][i].c << ' ';
            }
            //essentially transposing matrix to file write s.t.
            //vector downwards in space, with each column stepping in time
        }
        ufile << '\n';
        rhofile << '\n';
        pfile << '\n';
        efile << '\n';
        cfile << '\n';
        finalT << '\n';
    }
    ufile.close();
    rhofile.close();
    pfile.close();
    efile.close();
    cfile.close();
    finalT.close();
    double tol = 1e-3;
    if(fabs(u_max - u_min) < tol){ u_min += -0.1; u_max += 0.1;}
    if(fabs(rho_max - rho_min) < tol){ rho_min += -0.1; rho_max += 0.1;}
    if(fabs(p_max - p_min) < tol){ p_min += -0.1; p_max += 0.1;}
    if(fabs(e_max - e_min) < tol){ e_min += -0.1; e_max += 0.1;}
    if(fabs(c_max - c_min) < tol){ c_min += -0.1; c_max += 0.1;}    

    double arr[10] = {u_min, u_max, rho_min, rho_max, p_min, p_max, e_min, e_max, 
        c_min, c_max};
    output.assign(arr, &arr[10]);
    return output;
}

void writetoGNUplot(TestCase CASE, Vector domain, int num_snaps, Vector time_vec, 
        int I_0, int I_f){
    Vector space = CASE.space;
    string gpname = CASE.gpname;
    string datfile = CASE.datfile;
    string title = CASE.title;
    string giffile = CASE.giffile;
    int nT = CASE.nT;

    string u_datfile = "./data/u_" + datfile;
    string rho_datfile = "./data/rho_" + datfile;
    string p_datfile = "./data/p_" + datfile;
    string e_datfile = "./data/e_" + datfile;
    string c_datfile = "./data/c_" + datfile;

    ofstream file;
    file.open(gpname);
    double maxU = domain[1]; 
    double minU = domain[0];
    double maxX = *max_element(space.begin(), space.end());
    double minX = *min_element(space.begin(), space.end());
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal gif animate delay 8 \n";
    file << "n = " << nT << '\n';
    file << "set out \'" << giffile << "\' \n";
    file << "set key off \n";
    
    file << "set xrange[" << minX << ":" << maxX << "] \n";
    file << "do for [i=0:" << nT-1 << "] { \n    ";
    file << "j = i+2 \n";
    file << "set multiplot layout 2,3 \n";
    //rho plot:
    file << "set title \"" << "Density - rho" << "\" \n";
    file << "set yrange[" << domain[2] - 0.3*fabs(domain[3]-domain[2]) << \
        ":" << domain[3] + 0.3*fabs(domain[3]-domain[2]) << "] \n";
    file << "plot \"" << rho_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    //u plot:
    file << "set title \"" << "Velocity - u" << "\" \n";
    file << "set yrange[" << minU - 0.3*fabs(maxU-minU) << ":" << maxU + 0.3*fabs(maxU-minU) << "] \n";
    file << "plot \"" << u_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    //c plot:
    file << "set title \"" << "Concentration - c" << "\" \n";
    file << "set yrange[" << domain[8] - 0.3*fabs(domain[9]-domain[8]) << \
        ":" << domain[9] + 0.3*fabs(domain[9]-domain[8]) << "] \n";
    file << "plot \"" << c_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    //p plot:
    file << "set title \"" << "Pressure - p" << "\" \n";
    file << "set yrange[" << domain[4] - 0.3*fabs(domain[5]-domain[4]) << \
        ":" << domain[5] + 0.3*fabs(domain[5]-domain[4]) << "] \n";
    file << "plot \"" << p_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";
    //e plot:
    file << "set title \"" << "Internal energy - e" << "\" \n";
    file << "set yrange[" << domain[6] - 0.3*fabs(domain[7]-domain[6]) << \
        ":" << domain[7] + 0.3*fabs(domain[7]-domain[6]) << "] \n";
    file << "plot \"" << e_datfile << "\" using 1:j with linespoints pt 12 ps 0.5 \n";

    
    file << "unset multiplot \n } \n";
    file << "exit";
    file.close();

    //create time snapshots file
    ofstream snapfile;
    snapfile.open(title + "_snapshots.gp");
    snapfile << "#!/usr/local/bin/gnuplot -persist \n";
    snapfile << "set terminal svg \n";
    snapfile << "set xrange[0:1] \n"; 
    snapfile << "set out './visualisation/"<< CASE.title << " snapshot.svg\' \n";
    snapfile << "set key below samplen 2 spacing .7 font \",6\" \n";
    snapfile << "set multiplot layout 2,2 \n";
    snapfile << "set title \"Pressure - p\" \n";
    snapfile << "plot ";
    int index;
    for(int i = 1; i<num_snaps+1; ++i){
        index = I_0 + (int)(i*(I_f - I_0)/num_snaps);      //must multiply i*nT first due to 
                                                            //rounding down integerdivision
        snapfile << "'" << p_datfile << "\' using 1:" << index << " title 'T = " << \
           time_vec[index-1] << "\' with linespoints pt 6 ps 0.3 lc rgb \'blue\' ";
        if(i != num_snaps){
            snapfile << ", ";
        }else{
            snapfile << "\n";
        }
    }
    snapfile << "set title \"Velocity - u\" \n";
    snapfile << "plot ";
    for(int i = 1; i<num_snaps+1; ++i){
        index = I_0 + (int)(i*(I_f - I_0)/num_snaps); 
        snapfile << "'" << u_datfile << "\' using 1:" << index << " title 'T = " << \
           time_vec[index-1] << "\' with linespoints pt 6 ps 0.3 lc rgb \'blue\' ";
        if(i != num_snaps){
            snapfile << ", ";
        }else{
            snapfile << "\n";
        }
    }
    snapfile << "set title \"Density - rho\" \n";
    snapfile << "plot ";
    for(int i = 1; i<num_snaps+1; ++i){
        index = I_0 + (int)(i*(I_f - I_0)/num_snaps); 
        snapfile << "'" << rho_datfile << "\' using 1:" << index << " title 'T = " << \
           time_vec[index-1] << "\' with linespoints pt 6 ps 0.3 lc rgb \'blue\' ";
        if(i != num_snaps){
            snapfile << ", ";
        }else{
            snapfile << "\n";
        }
    }
    snapfile << "set title \"Concentration - c\" \n";
    //snapfile << "set arrow from 0.25,0.0 to 0.25,1.0 nohead lc rgb \'red\'\n";
    //snapfile << "set arrow from 0.75,0.0 to 0.75,1.0 nohead lc rgb \'red\'\n";
    snapfile << "plot ";
    for(int i = 1; i<num_snaps+1; ++i){
        index = I_0 + (int)(i*(I_f - I_0)/num_snaps); 
        snapfile << "'" << c_datfile << "\' using 1:" << index << " title 'T = " << \
           time_vec[index-1] << "\' with linespoints pt 6 ps 0.3 lc rgb \'blue\' ";
        if(i != num_snaps){
            snapfile << ", ";
        }else{
            snapfile << "\n";
        }
    }
    snapfile << "exit";

}

void GNUplot_finalT(TestCase CASE, string filetype){
    ofstream file;
    file.open(CASE.title + " finalT.gp");
    string datfile = "./data/finalT" + CASE.datfile;
    file << "#!/usr/local/bin/gnuplot -persist \n";
    file << "set terminal " << filetype << " \n";
    file << "set xrange[0:1] \n";
    file << "set out './visualisation/"<< CASE.title << " finalT." << filetype <<"\' \n";
    file << "set key off \n";
    file << "set multiplot layout 2,2 \n";
    file << "set title \"Density - rho\" \n";
    file << "plot '"<< datfile << "\' using 1:3 with points pt 6 ps 0.3 lc rgb \'blue\' \n";
    file << "set title \"Velocity - u\" \n";
    file << "plot '"<< datfile << "\' using 1:2 with points pt 6 ps 0.3 lc rgb \'blue\' \n";
    file << "set title \"Pressure - p\" \n";
    file << "plot '"<< datfile << "\' using 1:4 with points pt 6 ps 0.3 lc rgb \'blue\' \n"; 
    file << "set title \"Internal energy - e\" \n";
    file << "plot '"<< datfile << "\' using 1:5 with points pt 6 ps 0.3 lc rgb \'blue\' \n";     
    file << "exit";
}

//--------------------------------------------------------------(end) DATA AND VISUALISATION //

void tellmethings(TestCase CASE, bool MUSCL, string scheme, string BC, bool chemistry){
    cout << "|| " << CASE.title << " || " << endl;
    cout << "spatial discretisation = " << CASE.n - 2*CASE.nG << endl;
    cout << "simulated time = " << CASE.Tf << "s" << endl;
    cout << "HLLC solver- run as ";
    if(MUSCL){
        cout << "second order scheme: MUSCL, slope limiter = " << scheme << endl;
    }else{ cout << "first order scheme " << endl;
    }
    cout << "with " << BC << " boundary conditions" << endl;
    cout << "Chemistry is turned: ";
    if(chemistry){
        cout << "ON" << endl;
    }else{ cout << "OFF" << endl;
    }
    cout << "number of time steps in simulated time = " << CASE.nT << endl;
    cout << "5 data files generated in ./data/: \n  *" << CASE.datfile << "(x5 properties)\n";
    cout << "3 .gp files generated which can be compiled for visualisation: \n";
    cout << "for gif animation: \n  *" << CASE.gpname << endl;
    cout << "for final time solution image: \n  *" << CASE.title << " finalT.gp" << endl;
    cout << "for staggered snapshots image: \n  *" << CASE.title << "_snapshots.gp" << endl;
    cout << "compile with: \"chmod u+x <name.gp> \" " << \
        "then run \"./<name>.gp\", to generate file located in ./visualisation/" << endl;
    cout << "(compiling a gif with more than 300 time steps will take a super long time and is not recommended)\n";
    cout << "done." << endl;
}


int main(){
    
    //~~~~ USER DEFINED SETTINGS HERE ~~~~//
    TestCase CASE("Test1");            //"ShockTest1" is the studied full-system problem 
                                            //Toro tests and c-advection test also in shockTests.H file
    bool MUSCL = 1;                         //second order extension scheme: on/off
    string scheme = "SuperBee";             //slope limiter type - only "Superbee" implemented here
    string BC = "transmissive";             //Boundary conditions:
                                            //options:  "left reflective wall" or "transmissive"
    bool wall_correction = 0;               //wall heating correction: on/off 
    bool chemistry = 0;                     //chemistry switch: on/off
    int ns = 20;                            //number of chemistry sub-cycle steps
                                            //number of cells and final time set in shockTests.H
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    
    Vector time_vec;
    vector< vector<Prim> > MATRIX = generate_solution(CASE, time_vec, MUSCL, scheme, 
                                                BC, wall_correction, chemistry, ns);
    CASE.nT = MATRIX.size();

    tellmethings(CASE, MUSCL, scheme, BC, chemistry);  

    Vector domain = writetofile(CASE, MATRIX);
    
    //visualisation:
    //set parameters for desired snapshot image settings:
    int num_snaps = 10;     //number of time snap shots
    double Vt0 = 0.0;      //`visual t0' first snap shot time
    double Vtf = 1.0;      //'visual tf' final snap shot time
    //(ensure that these occur within the simulation run-time set in shockTests.H)

    int I_0, I_f;       //corresponding index in time vector
    //assert( Vt0 < Vtf && Vtf <= CASE.Tf);
    for(int i = 0; i < CASE.nT; ++i){
        if(time_vec[i] <= Vt0){
            I_0 = i;
        }
        if(time_vec[i] <= Vtf){
            I_f = i;
        }
    }

    writetoGNUplot(CASE, domain, num_snaps, time_vec, I_0, I_f);
    GNUplot_finalT(CASE, "svg");
 
    
    return 0;
}




