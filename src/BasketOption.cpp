//
// Created by MacBookPro on 12/09/2019.
//

#include <iostream>
#include "BasketOption.hpp"

BasketOption::BasketOption(double T_, int nbTimeSteps_, int size_, double strike_, PnlVect *weights_)
{
    this->T_ = T_;
    this->nbTimeSteps_ = nbTimeSteps_;
    this->size_ = size_;
    this->strike_ = strike_;
    this->weights_ = weights_;
    this->row = pnl_vect_create_from_zero(size_);
}

BasketOption::~BasketOption()
{
    pnl_vect_free(&weights_);
    pnl_vect_free(&row);
}

double BasketOption::payoff(const PnlMat *path)
{
    pnl_mat_get_row(row, path, nbTimeSteps_);
    double s_basket = pnl_vect_scalar_prod(row, weights_);
    return ((s_basket - strike_) > 0) ? (s_basket - strike_): 0;
}