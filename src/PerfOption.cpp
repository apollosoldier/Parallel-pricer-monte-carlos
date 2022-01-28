//
// Created by MacBookPro on 12/09/2019.
//

#include "PerfOption.hpp"
#include <iostream>

PerfOption::PerfOption(double T_, int nbTimeSteps_, int size_, PnlVect *weights)
{
    this->T_ = T_;
    this->nbTimeSteps_ = nbTimeSteps_;
    this->size_ = size_;
    this->weights_ = weights;
    this->row_i = pnl_vect_create_from_zero(size_);
    this->row_i_1 = pnl_vect_create_from_zero(size_);
}

PerfOption::~PerfOption()
{
    pnl_vect_free(&weights_);
    pnl_vect_free(&row_i);
    pnl_vect_free(&row_i_1);
}

double PerfOption::payoff(const PnlMat *path)
{
    double s_perf = 0;
    for (int i = 1 ; i<nbTimeSteps_+1; i++){
        pnl_mat_get_row(row_i, path, i);
        pnl_mat_get_row(row_i_1, path, i-1);
        double perf_denom  = pnl_vect_scalar_prod(row_i, weights_);
        double perf_nom = pnl_vect_scalar_prod(row_i_1, weights_);
        s_perf += (((perf_denom/perf_nom) - 1) > 0) ? ((perf_denom/perf_nom) - 1): 0;
    }
    return 1+s_perf;
}
