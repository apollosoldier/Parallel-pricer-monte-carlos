//
// Created by MacBookPro on 12/09/2019.
//

#ifndef MC_PRICER_PERFOPTION_HPP
#define MC_PRICER_PERFOPTION_HPP

#include <pnl/pnl_matvect.h>
#include "Option.hpp"

/// \brief Classe Option Performance

class PerfOption : public Option {
public:
    double strike_;
    PnlVect *weights_;
    PnlVect *row_i;
    PnlVect *row_i_1;

    virtual double payoff(const PnlMat *path);

    /**
     * Construit l'objet d'une Option Performance
     *
     * @param[in] T maturité de l'option
     * @param[in] nbTimeSteps_ nombre de time step
     * @param[in] size_ nombre des actifs de l'option
     * @param[in] Weights vecteur contenant les poids des actifs
     */
    PerfOption(double T_, int nbTimeSteps_, int size_, PnlVect *weights);

    virtual ~PerfOption();
};


#endif //MC_PRICER_PERFOPTION_HPP
