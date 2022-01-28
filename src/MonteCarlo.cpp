

#include "MonteCarlo.hpp"
#include <tgmath.h>
#include <iostream>
#include <mpi.h>
#include <limits>

#define MAX_SAMPLES 100000000;

MonteCarlo::MonteCarlo(BlackScholesModel *mod_, Option *opt_, PnlRng *rng_, double fdStep_, int nbSamples_) {
    this->mod_ = mod_;
    this->opt_ = opt_;
    this->rng_ = rng_;
    this->fdStep_ = fdStep_;
    this->nbSamples_ = nbSamples_;
    this->path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);//creating an empty matrix of size (nbTimeSteps+1) x d
}

void MonteCarlo::price(double &prix, double &ic) {
    // res is the estimated price of the option
    double res = 0;

    // simple variance of the estimator
    double var = 0;

    // generating nbSamples paths by BlackScholes model
    for (int i = 0; i < nbSamples_; i++) {

        //generating a path using BlackScholes
        this->mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);

        // computing the payOff option on this path
        res += opt_->payoff(path);
        var += pow(opt_->payoff(path), 2);
    }

    var = exp(-2 * mod_->r_ * opt_->T_) * ((var / nbSamples_) - pow(res / nbSamples_, 2));
    res *= exp(-mod_->r_ * opt_->T_) / nbSamples_;

    //returning the results
    prix = res;
    ic = 2 * 1.96 * sqrt(var) / sqrt(nbSamples_);
    pnl_mat_free(&path);
}


void MonteCarlo::price_parall_precision(int &hasFinished, int& opt_nbSamples, double &prix, double &ic, int mpi_size, int rank, double precision){
    double res=0, var=0, thread_res=0, thread_var=0;
    ic = 0;

    for (int i = 1; i < std::numeric_limits<double>::infinity(); i++) {

        //generating a path using BlackScholes
        this->mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);

        // computing the payOff option on this path
        thread_res += opt_->payoff(path);
        thread_var += pow(opt_->payoff(path), 2);
        MPI_Reduce(&thread_res, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&thread_var, &var, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (rank == 0){
            var = exp(-2 * mod_->r_ * opt_->T_) * ((var / (i*4)) - pow(res / (i * 4), 2));
            res *= exp(-mod_->r_ * opt_->T_) / (i*4);

            prix = res;
            ic = 2 * 1.96 * sqrt(var) / sqrt(i*4);
            if (ic/(2*1.96) < precision){
                opt_nbSamples = i * 4;
                hasFinished = 1;
                for (int j = 1; j < mpi_size; j++){
                    MPI_Send(&hasFinished, 1, MPI_INT, j,1, MPI_COMM_WORLD);
                }
                break;

            }
            else{
                for (int j = 1; j < mpi_size; j++){
                    MPI_Send(&hasFinished, 1, MPI_INT, j,1, MPI_COMM_WORLD);
                }
            }
        }
        else{
            MPI_Recv(&hasFinished, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
            if (hasFinished == 1)
                break;
        }


    }


}
void MonteCarlo::price_parall(double &prix, double &ic, int mpi_size, int rank){
    double res=0, var=0, thread_res=0, thread_var=0;

    price_slave((int) nbSamples_/mpi_size, thread_res, thread_var);
    //res += thread_res;
    MPI_Reduce(&thread_res, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //var += thread_var;
    MPI_Reduce(&thread_var, &var, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    //MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        price_master(res, var, prix, ic);
    }

}

void MonteCarlo::price_master(double &res, double &var, double &prix, double &ic){
    var = exp(-2 * mod_->r_ * opt_->T_) * ((var / nbSamples_) - pow(res / nbSamples_, 2));
    res *= exp(-mod_->r_ * opt_->T_) / nbSamples_;

    //returning the results
    prix = res;
    ic = 2 * 1.96 * sqrt(var) / sqrt(nbSamples_);
}

void MonteCarlo::price_slave(int samples_per_thread, double &thread_prix, double &thread_var){

    // generating nbSamples paths by BlackScholes model
    for (int i = 0; i < samples_per_thread; i++) {

        //generating a path using BlackScholes
        this->mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);

        // computing the payOff option on this path
        thread_prix += opt_->payoff(path);
        thread_var += pow(opt_->payoff(path), 2);
    }
}

void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic) {
    // res is the estimated price of the option
    double res = 0;

    // simple variance of the estimator
    double var = 0;

    //creating an empty matrix of size (nbTimeSteps+1) x d
    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    // generating nbSamples paths by BlackScholes model
    for (int i = 0; i < nbSamples_; i++) {
        //generating a path using BlackScholes
        this->mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);

        // computing the payOff option on this path
        res += opt_->payoff(path);
        var += pow(opt_->payoff(path), 2);

    }
    pnl_mat_free(&path);
    var = exp(-2 * mod_->r_ * (opt_->T_ - t)) * ((var / nbSamples_) - pow(res / nbSamples_, 2));
    res *= exp(-mod_->r_ * (opt_->T_ - t)) / nbSamples_;

    //returning the results
    prix = res;
    ic = 2 * 1.96 * sqrt(var) / sqrt(nbSamples_);
}

void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta, PnlVect *ic) {

    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
    PnlMat *upShifted = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);
    PnlMat *downShifted = pnl_mat_create(opt_->nbTimeSteps_ + 1, mod_->size_);

    // generating nbSamples paths by BlackScholes model
    for (int j = 0; j < nbSamples_; j++) {
        double res = 0;
        this->mod_->asset(path, t, opt_->T_, opt_->nbTimeSteps_, rng_, past);
        // Iterating over the assets
        for (int d = 0; d < mod_->size_; d++) {
            //Initializing the shifted paths
            pnl_mat_clone(upShifted, path);
            pnl_mat_clone(downShifted, path);

            //Generating the upShifted and the down Shifted paths
            this->mod_->shiftAsset(upShifted, path, d, fdStep_, t, opt_->T_ / opt_->nbTimeSteps_);
            this->mod_->shiftAsset(downShifted, path, d, -fdStep_, t, opt_->T_ / opt_->nbTimeSteps_);

            //Computing the difference

            double diffRes = (opt_->payoff(upShifted) - opt_->payoff(downShifted))/(2 * fdStep_ * MGET(past, past->m - 1, d));
            LET(delta, d) += (opt_->payoff(upShifted) - opt_->payoff(downShifted)) * exp(-mod_->r_ * (opt_->T_ - t)) / (nbSamples_ * 2 * fdStep_ * MGET(past, past->m - 1, d));
            LET(ic, d) += SQR(diffRes);
        }
    }
    pnl_vect_mult_scalar(ic, exp(-2*mod_->r_ * (opt_->T_ - t))/nbSamples_);
    for (int d = 0; d<mod_->size_; d++){
        LET(ic, d) -= SQR(GET(delta, d));
        LET(ic, d) = sqrt(GET(ic,d)/nbSamples_);
    }

    pnl_mat_free(&upShifted);
    pnl_mat_free(&downShifted);
    pnl_mat_free(&path);
}