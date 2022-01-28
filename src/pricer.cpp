
#include <iostream>
#include <string>
#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"
#include "VanillaCall.hpp"
#include "pnl/pnl_vector.h"
#include "jlparser/parser.hpp"
#include "AsianOption.hpp"
#include "PerfOption.hpp"
#include "BasketOption.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"



using namespace std;

int main(int argc, char **argv) {
    double T, r, strike, rho;
    PnlVect *spot, *sigma, *divid, *weights;
    string type;
    int size;
    int timeStep;
    size_t n_samples;

    char *infile = argv[1];
    Param *P = new Parser(infile);

    P->extract("option type", type);
    P->extract("maturity", T);
    P->extract("option size", size);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);
    if (type.compare("vanillacall")!=0)
        P->extract("payoff coefficients", weights, size);
    P->extract("interest rate", r);
    if (P->extract("dividend rate", divid, size, true) == false)
    {
        divid = pnl_vect_create_from_zero(size);
    }
    P->extract("correlation", rho);
    P->extract("timestep number", timeStep);
    if (type.compare("performance")!=0)
        P->extract("strike", strike);
    P->extract("sample number", n_samples);

    int nbSteps = timeStep;

    //Creating the option
    Option *op = nullptr;
    if (type.compare("asian")==0){
        op = new AsianOption(T,nbSteps, size, strike, weights);
    }
    else if (type.compare("basket")==0){
        op = new BasketOption(T,nbSteps, size, strike, weights);
    }
    else if (type.compare("performance")==0){
        op = new PerfOption(T, nbSteps, size, weights);
    }
    else{
        op = new VanillaCall(T, nbSteps, size, strike);
    }

    //Starting the parallelization *********************************
    clock_t debut, fin;
    debut = clock();

    int mpi_size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    //Initializing the random numbers generator
    PnlRng *rng_test = pnl_rng_dcmt_create_id(rank, time(0));
    //Initializing the generator's seed
    pnl_rng_sseed(rng_test, time(0));

    //Creating the BS model
    BlackScholesModel *bs_model = new BlackScholesModel(size,r,rho,sigma,spot, timeStep);

    //Creating the MonteCarlo simulator
    MonteCarlo  *mc = new MonteCarlo(bs_model, op, rng_test, 0.01, n_samples);
    double prix = 0;
    double ic = 0;

    //Version Parallele

    if (argv[2] != NULL)
    {
        int opt_nbSamples=0;
        double precision = atof(argv[2]);
        int hasFinished = 0;
        mc->price_parall_precision(hasFinished,opt_nbSamples, prix, ic, mpi_size, rank, precision);
        if (rank == 0) {
            cout << "Le nombre optimal de tirages pour la précision: " << precision << " est: " << opt_nbSamples
                 << endl;

        }
    } else {
        mc->price_parall(prix, ic, mpi_size, rank);
        if (rank == 0){
            cout <<"Le nombre de tirages: "<< n_samples << endl;

        }
    }

    //Finalizing the parallelization *********************************
    MPI_Finalize ();
    fin = clock();
    if (rank == 0)
    {
        cout <<"L'estimateur du prix en 0: "<< prix<<endl;
        cout << "La largeur de l'intervalle de confiance: " << ic/(2*1.96) << endl;
        cout << "Le temps de calcul: " << double(fin - debut) / CLOCKS_PER_SEC << endl;

    }

    //Version sequentielle
    //The option price:
    //mc->price(prix, ic);

    pnl_vect_free(&divid);
    delete bs_model;
    delete mc;
    delete op;
    pnl_rng_free(&rng_test);
    delete P;
}


