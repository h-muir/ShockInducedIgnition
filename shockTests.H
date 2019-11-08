//Module: Toro Test input cases

//non-dimensionalisation parameters:
//pre-shock values:
double rho1 = 1.2, u1 = 1024.7, p1 = 1e5;
//double rho2 = 4.6286, u2 = 265.667, p2 = 10.333e5;
double a0 = u1/3; //Mach 3 incident shock //t0=1.73016;
double L0 = a0*t0;
//shock ratios:
double R_rho = 0.259259, R_u = -1.35769, R_p = 0.0967742;
//post-shock values:
double rho2 = rho1/R_rho, u2 = 265.667, p2 = p1/R_p;

class TestCase{
    public:
        int n;          //spatial discretisation- number of points
        double L;       //1D spatial Length
        double Tf;      //final time
        double T0;      //initial time
        double CFL;     //Courant number
        double dx;      //spatial step
        double x0;      //property discontinuity position
        Prim WL;
        Prim WR;
        Vector space;
        string gpname, datfile, title, giffile;
        int nG;         //number of ghost cells at the boundary: on each side
                        //2 required for reflected shock BCs
        int nT;         //number of time steps,
                        //not assigned until after solution and number of time steps are known
        TestCase(string CASE):
            WL(0,0,0,0), 
            WR(0,0,0,0) {
                if(CASE == "ShockTest1"){
                    title = "Shock Test 1";
                    datfile = "shock_test1.dat";
                    gpname = "gifplot_shocktest1.gp";
                    giffile = "./visualisation/shocktest1.gif";
                    WL.rho = R_rho; WL.u = R_u; WL.p = R_p; WL.c = 1;
                    WR.rho = R_rho; WR.u = R_u; WR.p = R_p; WR.c = 1;
                    nG = 2;
                    n = 2000 + 2*nG;     //4 ghost cells for BCs: i = 0,1|---|n-2,n-1
                    L = 1;               //space domain 0 < x < 1
                    Tf = 1.40;           //simulation run to this final time
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.5;
                }

                if(CASE == "ShockTest2"){
                    title = "Shock Test 2";
                    datfile = "shock_test2.dat";
                    gpname = "gifplot_shocktest2.gp";
                    giffile = "./visualisation/shocktest2.gif";
                    WL.rho = 1.0; WL.u = 0.0; WL.p = 1.0; WL.c = 1;
                    WR.rho = R_rho; WR.u = R_u; WR.p = R_p; WR.c = 1;
                    nG = 2;
                    n = 100 + 2*nG;      //4 ghost cells for BCs: i = 0,1|---|n-2,n-1
                    L = 1;               //space domain 0 < x < 1
                    Tf = 1;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.03;
                }

                if(CASE == "TestC"){
                    title = "Test C";
                    datfile = "euler_eqs_testC.dat";
                    gpname = "gifplot_testC.gp";
                    giffile = "./visualisation/testC.gif";
                    WL.rho = 1.0; WL.u = 1.0; WL.p = 1.0; WL.c = 1.0;
                    WR.rho = 1.0; WR.u = 1.0; WR.p = 1.0; WR.c = 0.0;
                    nG = 2;
                    n = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.51;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.25;
                }

                if(CASE == "TestC2"){
                    title = "Test C2";
                    datfile = "euler_eqs_testC2.dat";
                    gpname = "gifplot_testC2.gp";
                    giffile = "./visualisation/testC2.gif";
                    WL.rho = 1.0; WL.u = 1.2; WL.p = 1.0; WL.c = 1.0;
                    WR.rho = 0.2; WR.u = 0.4; WR.p = 0.5; WR.c = 0.0;
                    nG = 2;
                    n = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.2;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.3;
                }
                if(CASE == "Test1"){
                    title = "Toro Test 1";
                    datfile = "euler_eqs_test1.dat";
                    gpname = "gifplot_test1.gp";
                    giffile = "./visualisation/torotest1.gif";
                    WL.rho = 1.0; WL.u = 0.75; WL.p = 1.0;
                    WR.rho = 0.125; WR.u = 0.0; WR.p = 0.1;
                    nG = 2;
                    n = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.2;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.3;
                }
                if(CASE == "Test2"){
                    title = "Toro412";
                    datfile = "euler_eqs_test2.dat";
                    gpname = "gifplot_test2.gp";
                    giffile = "./visualisation/torotest2.gif";
                    WL.rho = 1.0; WL.u = -2.0; WL.p = 0.4;
                    WR.rho = 1.0; WR.u = 2.0; WR.p = 0.4;
                    nG = 2;
                    n = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.15;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.5;
                }
                if(CASE == "Test3"){
                    title = "Toro413";
                    datfile = "euler_eqs_test3.dat";
                    gpname = "gifplot_test3.gp";
                    giffile = "./visualisation/torotest3.gif";
                    WL.rho = 1.0; WL.u = 0.0; WL.p = 1000.0;
                    WR.rho = 1.0; WR.u = 0.0; WR.p = 0.01;
                    nG = 2;
                    n = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.012;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.5;
                }
                if(CASE == "Test4"){
                    title = "Toro Test 4";
                    datfile = "euler_eqs_test4.dat";
                    gpname = "gifplot_test4.gp";
                    giffile = "./visualisation/torotest4.gif";
                    WL.rho = 5.99924; WL.u = 19.5975; WL.p = 460.894;
                    WR.rho = 5.99242; WR.u = -6.19633; WR.p = 46.0950;
                    nG = 2;
                    n = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.035;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.4;
                }
                if(CASE == "Test5"){
                    title = "Toro Test 5";
                    datfile = "euler_eqs_test5.dat";
                    gpname = "gifplot_test5.gp";
                    giffile = "./visualisation/torotest5.gif";
                    WL.rho = 1.0; WL.u = -19.59745; WL.p = 1000.00;
                    WR.rho = 1.0; WR.u = -19.59745; WR.p = 0.01;
                    nG = 2;
                    n = 100 + 2*nG;
                    L = 1.0;    //space domain 0 < x < 1
                    Tf = 0.012;
                    T0 = 0;
                    CFL = 0.9;
                    dx = L/(n-1-2*nG);
                    x0 = 0.8;
                }

                space.resize(n);
                for(int i = 0; i<n; ++i){
                    space[i] = i*dx - nG*dx;
                }
                
            }
};


