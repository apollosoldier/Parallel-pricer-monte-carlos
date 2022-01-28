//
// Created by MacBookPro on 12/09/2019.
//

#ifndef MC_PRICER_BASKETOPTION_HPP
#define MC_PRICER_BASKETOPTION_HPP


#include <pnl/pnl_matvect.h>
#include "Option.hpp"

/// \brief Classe Option Basket


class BasketOption : public Option{
public:
    double strike_;
    PnlVect *weights_;
    PnlVect *row;

    virtual double payoff(const PnlMat *path);

    /**
     * Construit l'objet d'une Option Performance
     *
     * @param[in] T maturité de l'option
     * @param[in] nbTimeSteps_ nombre de time step
     * @param[in] size_ nombre des actifs de l'option
     * @param[in] strike_ prix de l'exercice de l'option
     * @param[in] Weights vecteur contenant les poids des actifs
     */
     BasketOption(double T_, int nbTimeSteps_, int size_, double strike_, PnlVect *weights);

     virtual ~BasketOption();
};


#endif //MC_PRICER_BASKETOPTION_HPP
