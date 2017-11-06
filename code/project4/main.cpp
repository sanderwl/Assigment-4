#include <iostream>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> uniformDist(0.0,1.0);

ofstream outFile, outFile2;

double random(){
    return uniformDist(randomEngine);
}

void toFile(int cycles, double Etemp){
    outFile << cycles << ", " << Etemp << endl;
}

mat randomMatrix(mat &A, int L){

    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            if(random() < 0.5){
                A(i,j) = 1;
            }
            else{
                A(i,j) = -1;
            }
        }
    }
    return A;
}


void MP(int L, mat &A, vec &prob, int &acceptance, double &E, double &Mtemp, double &E_2, double &Etemp, int cycle, double &M, double &M_2){
    int xp, xn, yp, yn;
    double deltaE;
    for(int x = 0; x < L; x++){
        for(int y = 0; y < L; y++){

            xp = (x + 1) % L;
            yp = (y + 1) % L;
            xn = (x - 1 + L) % L;
            yn = (y - 1 + L) % L;

            deltaE = 2*A(x,y) * (  A(xp,y) + A(xn,y) + A(x, yn) + A(x, yp));

            if(random() <= prob(deltaE + 8)){
                A(x, y) *= -1;
                Etemp += deltaE;

                Mtemp += 2*A(x,y); //Mean magnetic moment
                acceptance += 1;

            }
        }
    }
    E += Etemp;
    E_2 += Etemp*Etemp;
    M += (Mtemp);
    M_2 += Mtemp*Mtemp;
    toFile(cycle, Etemp);

    //cout << Etemp << endl;
}


void openFiles(){
    string outFileName = "mc_cycles.txt";
    string outFileName2 = "averegeE.txt";
    outFile.open(outFileName);
    outFile2.open(outFileName2);
}

int main()
{

    double k = 1.0; //J/K
    double T = 1.0; //kT/J
    double beta = 1.0/(k*T);
    double J = 1.0;

    int mcs = 100000;

    int L = 2;

    openFiles();

    //mat A = zeros(L,L);
    //randomMatrix(A, L);
    //cout << A << endl;
    mat A = ones(L,L);

    double Etemp = 0;

    //Initial values of the temporary energy
    for(int x = 0; x<L; x++){
        for(int y = 0; y<L;y++){

            int xn = (x - 1 + L) % L;
            int yn = (y - 1 + L) % L;

            Etemp -= A(x,y) * (  A(xn,y) + A(x, yn));
        }
    }

    double E = 0;
    double E_2 = 0;
    double Mtemp = 0;
    double M = 0;
    double M_2= 0;

    vec prob(17);
    for(int i=-8; i <= 8; i+=4){
        prob(i+8) = 0;
    }
    for(int i=-8; i <= 8; i+=4){
        prob(i + 8) = exp(-i/T);
    }
    //cout << prob << endl;
    int acceptance = 0;
    for(int cycles=0;cycles<=mcs;cycles++){
        MP(L, A, prob, acceptance, E, Mtemp, E_2, Etemp, cycles, M, M_2);
        toFile(cycles, Etemp);

    }
    double averegeE = E/(L*L*mcs); //Average energy
    double averegeESquared = E_2/(L*L*mcs);

    //double magnetic = accu(A); //Magnetic moment
    double averegeM = (M)/(mcs);
    double averegeMSquared = (M_2/(mcs));

    double PF = 2*exp(-8*J*beta) + 2*exp(8*J*beta) + 12; //Partition function
    //double sus = beta*((M_2)/(PF) - (M*M)/(PF*PF));
    double sus = beta*((averegeMSquared) - ((averegeM))*((averegeM)));//(L*L);

    cout << averegeE<< endl;
    cout << averegeM << " " << averegeMSquared << endl;
    cout << sus << endl;


    outFile.close();
    outFile2.close();

    return 0;
}
