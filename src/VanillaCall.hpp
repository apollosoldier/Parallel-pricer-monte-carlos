//
// Created by MacBookPro on 12/09/2019.
//

#ifndef MC_PRICER_VANILLACALL_HPP
#define MC_PRICER_VANILLACALL_HPP

#include "Option.hpp"

/// \brief Classe Option Call Vanille

class VanillaCall: public Option {
    public:
        double strike_;

        VanillaCall(double T_, int nbTimeSteps_, int size_, double strike_);

        virtual double payoff(const PnlMat *path);
        virtual ~VanillaCall(){};
};


#endif //MC_PRICER_VANILLACALL_HPP
