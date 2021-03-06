#include <iostream>
#include <armadillo>
#include <random>

using namespace std;
using namespace arma;

random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> uniformDist(0.0,1.0);

ofstream outFile, outFile2, outFile3, outFile4;

//Using the c++ RNG
double random(){
    return uniformDist(randomEngine);
}
//Printing magnetization, energy and acceptance to file
void toFile(double M, double E, double T, int acceptance){
    outFile << M << ", " << endl;
    outFile2 << E << ", " << endl;
    outFile4 << acceptance << ", " << endl;
}
//Printing total number of acceptance per temperature value
void toFile2(int acceptance, double T){
    outFile3 << T << ", " << acceptance << endl;
}
//Creating a random matrix
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

//Metropolis algorithm
void MP(int L, mat &A, vec prob, int &acceptance, double &E, double &Mtemp, double &E_2, double &Etemp, double &M, double &M_2){
    int xp, xn, yp, yn;
    double deltaE;
    for(int x = 0; x < L; x++){
        for(int y = 0; y < L; y++){

            int xrand = round(random()*(L-1));
            int yrand = round(random()*(L-1));

            xp = (xrand + 1) % L;
            yp = (yrand + 1) % L;
            xn = (xrand - 1 + L) % L;
            yn = (yrand - 1 + L) % L;

            deltaE = 2.0 * A(xrand,yrand) * (  A(xp,yrand) + A(xn,yrand) + A(xrand, yn) + A(xrand, yp));

            if(random() <= prob(deltaE + 8)){
                A(xrand, yrand) *= -1;
                Etemp += deltaE;

                Mtemp += 2*A(xrand,yrand);
                acceptance += 1;

            }
        }
    }
    E += Etemp;
    E_2 += Etemp*Etemp;
    M += abs(Mtemp);
    M_2 += Mtemp*Mtemp;
}

//Opening files
void openFiles(double T){
        string outFileName = "magnetization" + to_string(T) + ".txt";
        string outFileName2 = "energy" + to_string(T) + ".txt";
        string outFileName4 = "acceptanceMC" + to_string(T) + ".txt";
        outFile.open(outFileName);
        outFile2.open(outFileName2);
        outFile4.open(outFileName4);
}

void openFiles2(){
    string outFileName3 = "TempAcceptance.txt";
    outFile3.open(outFileName3);
}

int main()
{
    //openFiles2();
    for(double T = 1.0; T <= 2.4; T+=1.4){
    double k = 1.0; //J/K
    double beta = 1.0/(k*T);

    int mcs = 100000; //Total number of Monte Carlo cycles

    int L = 2; //Lattice size

    //openFiles(T);

    mat A = ones(L,L);
    randomMatrix(A, L); //Comment out to have an initial matrix with only up spins
    //cout << A << endl;

    //Initial values of the temporary energy
    double Etemp = 0;
    for(int x = 0; x<L; x++){
        for(int y = 0; y<L;y++){

            int xn = (x - 1 + L) % L;
            int yn = (y - 1 + L) % L;

            Etemp -= A(x,y) * (  A(xn,y) + A(x, yn)); //Calculating the inital energy value
        }
    }

    double E = 0;
    double E_2 = 0;
    double Mtemp = accu(A); //Initial magnetization value
    double M = 0;
    double M_2= 0;

    //Boltzmann probability
    vec prob(17);
    for(int i=-8; i <= 8; i+=4){
        prob(i+8) = 0;
    }
    for(int i=-8; i <= 8; i+=4){
        prob(i + 8) = exp(-i/T);
    }

    int acceptance = 0;
    //Monte Carlo loop
    for(int cycles=0;cycles<=mcs;cycles++){
        MP(L, A, prob, acceptance, E, Mtemp, E_2, Etemp, M, M_2);


        //if(cycles >= 100000){ //Only values after equilibrum
            //int n=1;
            //if((cycles % n) == 0){ //Only writing every n-th output to file
                //toFile(Mtemp, Etemp, T, acceptance); //Divide by "cycles" to get the mean M and mean E
            //}
        //}
    }

    //toFile2(acceptance, T);

    double averegeE = E/(mcs); //Average energy
    double averegeESquared = E_2/(mcs); //Average energy squared

    double averegeM = (M)/(mcs); //Average magnetization
    double averegeMSquared = (M_2/mcs); //Average magnetization squared

    //double PF = 2*exp(-8*J*beta) + 2*exp(8*J*beta) + 12; //Partition function
    double sus = beta*(averegeMSquared - averegeM*averegeM); //Susceptibility

    double heat = (beta*(averegeESquared - averegeE*averegeE))/T; //Specific heat

    cout << endl << "Average energy: " << averegeE << " while T = " << to_string(T) << endl;
    cout << "Average magnetization: " << averegeM << " while T = " << to_string(T) << endl;
    cout << "Specific heat: " << heat << " while T = " << to_string(T) << endl;
    cout << "Susceptibility: " << sus << " while T = " << to_string(T) << endl;
    cout << "Variance: " << averegeESquared - averegeE*averegeE << " while T = " << to_string(T) << endl;

    //outFile.close();
    //outFile2.close();
    //outFile4.close();
    }
    //outFile3.close();



    return 0;
}
