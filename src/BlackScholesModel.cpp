
#include "BlackScholesModel.hpp"
#include <math.h>
#include <ctime>
#include <iostream>
#include <assert.h>

using namespace std;

BlackScholesModel::BlackScholesModel()
{
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_vect_free(&sigma_);
    pnl_vect_free(&spot_);
    pnl_mat_free(&correlationMatrix);
    pnl_mat_free(&G);
    pnl_vect_free(&g);
    pnl_vect_free(&lineD);
}

BlackScholesModel::BlackScholesModel(int size_, double r_, double rho_, PnlVect *sigma_, PnlVect *spot_, PnlVect *trend_, int nbTimeSteps):BlackScholesModel(size_, r_, rho_, sigma_, spot_, nbTimeSteps)
{
    this->trend_ = trend_;
}

BlackScholesModel::BlackScholesModel(int size_, double r_, double rho_, PnlVect *sigma_, PnlVect *spot_, int nbSteps)
{
    this->size_ = size_;
    this->r_ = r_;
    this->rho_ = rho_;
    this->sigma_ = sigma_;
    this->spot_ = spot_;
    this->correlationMatrix = pnl_mat_create_from_scalar(size_, size_, rho_);

    //Setting the diag of the correlationMatrix to 1
    for (int i = 0; i < size_; i++){
        MLET(correlationMatrix, i,i) = 1;
    }

    // Computing Cholesky decomposition of the correlationMatrix: L
    int cholRes = pnl_mat_chol(correlationMatrix);

    //Creating the matrix G
    this->G = pnl_mat_create_from_zero(size_,nbSteps+1);

    //Creating line D
    lineD = pnl_vect_create_from_zero(size_);

    //Creating vect g
    g = pnl_vect_create_from_zero(size_);

    PnlVect* v;
}


void BlackScholesModel::asset(PnlMat *path, double t, double T, int nbTimeSteps, PnlRng *rng, const PnlMat *past)
{
    // If pricing is done at time 0
    if(t==0)
    {
        this->asset(path, T, nbTimeSteps, rng);
        return;
    }

    assert(past->m <= path->m);

    //Generating the Matrix G
    pnl_mat_rng_normal(G, size_, nbTimeSteps+1, rng);

    double step = T/nbTimeSteps;

    //lastTime is the last discretisation time in past
    int k = past->m - 2;
    double lastTime = k * step;

    for (int i = 0; i < nbTimeSteps + 1 ; i += 1){
        pnl_mat_get_col(g, G, i);
        for (int d = 0; d < size_; d++){
            if( i <= k ){
                MLET(path, i, d) = MGET(past, i, d);
            }
            else if ( i == k+1 && t - lastTime != step){
                // evaluation of the Stk+1 different than the rest
                double assetPrice = nextPrice(r_, MGET(past, k+1, d), lastTime + step - t, d, g);
                MLET(path, i, d) = assetPrice;
            }
            else{
                double assetPrice = nextPrice(r_, MGET(path, i-1, d), step, d, g);
                MLET(path, i, d) = assetPrice;
            }
        }
    }
}


// function created to factorize the code : compute the next price knowing the last one and the step
double BlackScholesModel::nextPrice(double trendOrRate, double lastPrice, double step, int d, PnlVect *g)
{
    //Extraction of the line number d from the decomposed correlation matrix: L
    pnl_mat_get_row(lineD, correlationMatrix, d);

    //Computing the asset price
    double term1 = (trendOrRate- SQR(GET(sigma_, d)) / 2) * step;
    double term2 = GET(sigma_,d) * sqrt(step) * pnl_vect_scalar_prod(lineD, g);
    double res = lastPrice * exp(term1 + term2);
    return res;
}


void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng)
{
    //Generate Matrix G
    pnl_mat_rng_normal(G, size_, nbTimeSteps+1, rng);

    //Computing the time step
    double step = T/nbTimeSteps;

    for (int d = 0; d < size_; d++){
        MLET(path, 0, d) = GET(spot_,d);
    }

    for (int t = 1; t < nbTimeSteps + 1 ; t += 1){
        pnl_mat_get_col(g, G, t);
        for (int d = 0; d < size_; d++){
            //Extraction of the line number d from the decomposed correlation matrix: L
            pnl_mat_get_row(lineD, correlationMatrix, d);

            //Computing the asset price
            double term1 = (r_- SQR(GET(sigma_, d)) / 2) * step;
            double term2 = GET(sigma_,d) * sqrt(step) * pnl_vect_scalar_prod(lineD, g);
            double assetPrice = MGET(path, t-1, d) * exp(term1 + term2);

            //Storing the asset price in path.
            MLET(path, t, d) = assetPrice;
        }
    }
}

void BlackScholesModel::shiftAsset(PnlMat *shift_path, const PnlMat *path, int d, double h, double t, double timestep)
{
    int index = ((floor(t / timestep))==ceil(t / timestep) )?  t/timestep : floor(t/timestep) + 1;
    for (int i = index ; i<path->m; i++){
        MLET(shift_path, i, d) = (1+h)*MGET(path,i, d);
    }
}



void BlackScholesModel::simul_market(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng)
{
//Generate Matrix G
    pnl_mat_rng_normal(G, size_, nbTimeSteps+1, rng);

    //Computing the time step
    double step = T/nbTimeSteps;

    for (int d = 0; d < size_; d++){
        MLET(path, 0, d) = GET(spot_,d);
    }

    for (int t = 1; t < nbTimeSteps + 1 ; t += 1){
        pnl_mat_get_col(g, G, t);
        for (int d = 0; d < size_; d++){
            //Storing the asset price in path.
            MLET(path, t, d) = nextPrice(GET(trend_, d), MGET(path, t-1, d), step, d,g);
        }
    }
}
