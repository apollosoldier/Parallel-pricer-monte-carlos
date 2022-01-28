//
// Created by MacBookPro on 12/09/2019.
//

#include "AsianOption.hpp"
#include <iostream>

AsianOption::AsianOption(double T_, int nbTimeSteps_, int size_, double strike_, PnlVect *weights)
{
    this->T_ = T_;
    this->nbTimeSteps_ = nbTimeSteps_;
    this->size_ = size_;
    this->strike_ = strike_;
    this->weights_ = weights;
    this->col_d = pnl_vect_create_from_zero(nbTimeSteps_+1);
}

AsianOption::~AsianOption()
{
    pnl_vect_free(&weights_);
    pnl_vect_free(&col_d);
}

double AsianOption::payoff(const PnlMat *path)
{
    double s_asian = 0;
    for(int d = 0; d<size_; d++){
        pnl_mat_get_col(col_d, path, d);
        s_asian +=GET(weights_, d) * pnl_vect_sum(col_d);
    }
    s_asian /= (nbTimeSteps_+1);
    return ((s_asian - strike_) > 0) ? (s_asian - strike_): 0;
}
