#include <iostream>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> uniformDist(0.0,1.0);

ofstream outFile, outFile2, outFile3, outFile4, outFile5, outFile6;

double random(){
    return uniformDist(randomEngine);
}

void toFile(double Mtemp, double Etemp, double T, int acceptance){
    if(T == 1.0){
        outFile << Mtemp << ", " << endl;
        outFile2 << Etemp << ", " << endl;
        outFile5 << acceptance << ", " << endl;
    }
    else if(T == 2.4){
        outFile3 << Mtemp << ", " << endl;
        outFile4 << Etemp << ", " << endl;
        outFile5 << acceptance << ", " << endl;
    }

}

mat randomMatrix(mat &A, int L){

    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            if(random() <= 0.5){
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
}


void openFiles(double T){
    if(T == 1.0){
        string outFileName = "magnetizationT1.txt";
        string outFileName2 = "energyT1.txt";
        string outFileName5 = "acceptanceT1.txt";
        outFile.open(outFileName);
        outFile2.open(outFileName2);
        outFile5.open(outFileName5);

    }
    else if(T == 2.4){
        string outFileName3 = "magnetizationT24.txt";
        string outFileName4 = "energyT24.txt";
        string outFileName6 = "acceptanceT24.txt";
        outFile3.open(outFileName3);
        outFile4.open(outFileName4);
        outFile6.open(outFileName6);
    }
}

int main()
{

    double k = 1.0; //J/K
    double T = 2.4; //kT/J
    double beta = 1.0/(k*T);
    double J = 1.0;

    int mcs = 1000;

    int L = 20;

    openFiles(T);

    mat A = ones(L,L);
    randomMatrix(A, L);
    //cout << A << endl;

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

        int n=100;
        if((cycles % n) == 0){
            toFile(Mtemp, Etemp, T, acceptance);
        }

    }
    double averegeE = E/(mcs); //Average energy
    double averegeESquared = E_2/(mcs);

    double averegeM = (M)/(mcs);
    double averegeMSquared = (M_2/mcs);

    double PF = 2*exp(-8*J*beta) + 2*exp(8*J*beta) + 12; //Partition function
    double sus = beta*(averegeMSquared - averegeM*averegeM);

    double heat = (beta*(averegeESquared - (averegeE*averegeE)))/T;

    cout << averegeE << " " << averegeESquared << endl;
    cout << averegeM << endl;
    cout << heat << endl;
    cout << sus << endl;

    outFile.close();
    outFile2.close();
    outFile3.close();
    outFile4.close();
    outFile5.close();
    outFile6.close();

    return 0;
}
