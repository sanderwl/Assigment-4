#include <iostream>
#include <armadillo>
#include <random>
#include <openmpi-x86_64/mpi.h>
#include "time.h"

using namespace std;
using namespace arma;

random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> uniformDist(0.0,1.0);

ofstream outFile, outFile2, outFile3, outFile4, outFile5;

int numProcs, myRank;
clock_t start,finish;

double rrandom(){
    return uniformDist(randomEngine);
}

void toFile(double Mtemp, double Etemp, double T, int acceptance){
    outFile << Mtemp << ", " << endl;
    outFile2 << Etemp << ", " << endl;
    outFile4 << acceptance << ", " << endl;
}

void toFile2(int acceptance, double T){
    outFile3 << T << ", " << acceptance << endl;
}


mat randomMatrix(mat &A, int L){

    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            if(rrandom() <= 0.5){
                A(i,j) = 1;
            }
            else{
                A(i,j) = -1;
            }
        }
    }
    return A;
}


void MP(int L, mat &A, vec prob, int &acceptance, double &E, double &Mtemp, double &E_2, double &Etemp, double &M, double &M_2){
    int xp, xn, yp, yn;
    double deltaE;
    for(int x = 0; x < L; x++){
        for(int y = 0; y < L; y++){

            int xrand = round(rrandom()*(L-1));
            int yrand = round(rrandom()*(L-1));

            xp = (xrand + 1) % L;
            yp = (yrand + 1) % L;
            xn = (xrand - 1 + L) % L;
            yn = (yrand - 1 + L) % L;

            deltaE = 2.0 * A(xrand,yrand) * (  A(xp,yrand) + A(xn,yrand) + A(xrand, yn) + A(xrand, yp));

            if(rrandom() <= prob(deltaE + 8)){
                A(xrand, yrand) *= -1;
                Etemp += deltaE;

                Mtemp += 2*A(xrand,yrand); //Mean magnetic moment
                acceptance += 1;

            }
        }
    }

}


void openFiles(){
        string outFileName = "energy100.txt";
        string outFileName2 = "magnetic100.txt";
        string outFileName3 = "heat100.txt";
        string outFileName4 = "sus100.txt";
        string outFileName5 = "temp.txt";
        outFile.open(outFileName);
        outFile2.open(outFileName2);
        outFile3.open(outFileName3);
        outFile4.open(outFileName4);
        outFile5.open(outFileName5);
}
void toFileBig(double averegeE, double averegeM, double heat, double sus, int L, int mcs, double T){

    outFile << averegeE << ", " << endl;
    outFile2 << averegeM << ", " << endl;
    outFile3 << heat << ", " << endl;
    outFile4 << sus << ", " << endl;
    outFile5 << T  << ", " << endl;

}

void openFiles2(){
    string outFileName3 = "TempAcceptance.txt";
    outFile3.open(outFileName3);
}

int main(int argn, char* argv[])
{
    MPI_Init(&argn, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&numProcs);

    start = clock();
    if(myRank == 0){
        openFiles();
    }
    //openFiles2();
    for(double T = 2.0; T <= 3.0; T+=0.05){
    double k = 1.0; //J/K
    double beta = 1.0/(k*T);

    int mcs = 1000000;

    int L = 100;


    mat A = ones(L,L);
    randomMatrix(A, L);
    //cout << A << endl;

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

//        int n=1;
//        if((cycles % n) == 0){
//            toFile(Mtemp, Etemp, T, acceptance);
//        }
        E += Etemp;
        E_2 += Etemp*Etemp;
        M += abs(Mtemp);
        M_2 += Mtemp*Mtemp;

    }

    MPI_Allreduce(&E,&E,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&E_2,&E_2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&M,&M,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&M_2,&M_2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    //toFile2(acceptance, T);

    double averegeE = E/(mcs*numProcs); //Average energy
    double averegeESquared = E_2/(mcs*numProcs);

    double averegeM = M/(mcs*numProcs);
    double averegeMSquared = M_2/(mcs*numProcs);

    //double PF = 2*exp(-8*J*beta) + 2*exp(8*J*beta) + 12; //Partition function
    double sus = beta*(averegeMSquared - averegeM*averegeM);

    double heat = (beta*(averegeESquared - averegeE*averegeE))/T;

    if(myRank == 0) {
        toFileBig(averegeE, averegeM, heat, sus, L, mcs, T);
    }

/*
    cout << endl << "Average energy: " << averegeE << " while T = " << to_string(T) << endl;
    cout << "Average magnetization: " << averegeM << " while T = " << to_string(T) << endl;
    cout << "Specific heat: " << heat << " while T = " << to_string(T) << endl;
    cout << "Susceptibility: " << sus << " while T = " << to_string(T) << endl;
    cout << "Variance: " << averegeESquared - averegeE*averegeE << " while T = " << to_string(T) << endl;*/
    if(myRank == 0){
        cout << "Temperature = " << to_string(T) << " and total T = 2.3" << endl;
        cout << sus << " " << heat << " " << averegeE << " " << averegeM << endl;
    }
    }

    if(myRank == 0) {
        outFile.close();
        outFile2.close();
        outFile3.close();
        outFile4.close();
        outFile5.close();
    }
    /*
    if(myRank == 0){
        outFile7.close();
    }
    */
    finish = clock();
    double time = (double (finish)- double (start))/(CLOCKS_PER_SEC);
    cout << time << endl;
    MPI_Finalize();

    return 0;
}
