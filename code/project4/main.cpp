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

void toFile(double Etemp){
    outFile2 << Etemp << ", " << endl;
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


void MP(int L, mat &A, vec &prob, int &acceptance, double &E, double &Mtemp, double &E_2, double &Etemp, double &M, double &M_2){
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
    M += abs(Mtemp);
    M_2 += Mtemp*Mtemp;
    //toFile(cycle, Etemp);
}


void openFiles(){
    string outFileName = "mc_cycles.txt";
    string outFileName2 = "energy.txt";
    outFile.open(outFileName);
    outFile2.open(outFileName2);
}

int main()
{

    double k = 1.0; //J/K
    double T = 1.0; //kT/J
    double beta = 1.0/(k*T);
    double J = 1.0;

    int mcs = 1000000;

    int L = 20;

    openFiles();

    mat A = zeros(L,L);
    randomMatrix(A, L);
    //mat A = ones(L,L);

    //Initial values of the temporary energy
    double Etemp = 0;
    for(int x = 0; x<L; x++){
        for(int y = 0; y<L;y++){

            int xn = (x - 1 + L) % L;
            int yn = (y - 1 + L) % L;

            Etemp -= A(x,y) * (  A(xn,y) + A(x, yn));
        }
    }

    double E = 0;
    double E_2 = 0;
    double Mtemp = accu(A);
    double M = 0;
    double M_2= 0;

    vec prob(17);
    for(int i=-8; i <= 8; i+=4){
        prob(i+8) = 0;
    }
    for(int i=-8; i <= 8; i+=4){
        prob(i + 8) = exp(-i/T);
    }

    int acceptance = 0;
    for(int cycles=0;cycles<=mcs;cycles++){
        MP(L, A, prob, acceptance, E, Mtemp, E_2, Etemp, M, M_2);
        toFile(Etemp);

    }
    double averegeE = E/(mcs); //Average energy
    double averegeESquared = E_2/(mcs);

    double averegeM = (M)/(mcs);
    double averegeMSquared = (M_2/mcs);

    double PF = 2*exp(-8*J*beta) + 2*exp(8*J*beta) + 12; //Partition function
    double sus = beta*(averegeMSquared - averegeM*averegeM);

    double heat = (beta/T)*(averegeESquared - (averegeE*averegeE));

    cout << averegeE << " " << averegeESquared << endl;
    cout << averegeM << endl;
    cout << heat << endl;
    cout << sus << endl;


    outFile.close();
    outFile2.close();

    return 0;
}
