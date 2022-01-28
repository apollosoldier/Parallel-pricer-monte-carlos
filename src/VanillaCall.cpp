//
// Created by MacBookPro on 12/09/2019.
//

#include "VanillaCall.hpp"

VanillaCall::VanillaCall(double T_, int nbTimeSteps_, int size_, double strike_)
{
    this->T_ = T_;
    this->nbTimeSteps_ = nbTimeSteps_;
    this->size_ = size_;
    this->strike_ = strike_;
}

double VanillaCall::payoff(const PnlMat *path)
{
    return ((MGET(path, nbTimeSteps_, 0) - strike_) > 0) ? (MGET(path, nbTimeSteps_, 0) - strike_): 0;
}