#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Modèle de Black Scholes
class BlackScholesModel
{
public:
    int size_; /// nombre d'actifs du modèle
    double r_; /// taux d'intérêt
    double rho_; /// paramètre de corrélation
    PnlVect *sigma_; /// vecteur de volatilités
    PnlVect *spot_; /// valeurs initiales des sous-jacents
    PnlMat *correlationMatrix;
    PnlMat *G;
    PnlVect * lineD;
    PnlVect * g;
    PnlVect *trend_;

    BlackScholesModel();

    /**
     * Constructeur par arguments
     *
     * @param[in] size_ La taille du modèle
     * @param[in] r_ Le taux d'interet
     * @param[in] rho_ Le coefficient de correlation
     * @param[in] sigma_ Le vecteur des volatilités
     * @param[in] spot_ Les prix initials
     * @param[in] nbTimeSteps Le nombres des dates de constatation
     */
    BlackScholesModel(int size_, double r_, double rho_, PnlVect *sigma_, PnlVect *spot_, int nbTimeSteps);

    /**
     * Constructeur par arguments
     *
     * @param[in] size_ La taille du modèle
     * @param[in] r_ Le taux d'interet
     * @param[in] rho_ Le coefficient de correlation
     * @param[in] sigma_ Le vecteur des volatilités
     * @param[in] spot_ Les prix initials
     * @param[in] trend_ vecteurs des tendances du marché
     * @param[in] nbTimeSteps Le nombres des dates de constatation
     */
    BlackScholesModel(int size_, double r_, double rho_, PnlVect *sigma_, PnlVect *spot_, PnlVect *trend_, int nbTimeSteps);

    ~BlackScholesModel();

    /**
     * Génère une trajectoire du modèle et la stocke dans path
     *
     * @param[out] path contient une trajectoire du modèle.
     * C'est une matrice de taille (nbTimeSteps+1) x d
     * @param[in] T  maturité
     * @param[in] nbTimeSteps nombre de dates de constatation
     */
    void asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);

    /**
     * Calcule une trajectoire du modèle connaissant le
     * passé jusqu' à la date t
     *
     * @param[out] path  contient une trajectoire du sous-jacent
     * donnée jusqu'à l'instant t par la matrice past
     * @param[in] t date jusqu'à laquelle on connait la trajectoire.
     * t n'est pas forcément une date de discrétisation
     * @param[in] nbTimeSteps nombre de pas de constatation
     * @param[in] T date jusqu'à laquelle on simule la trajectoire
     * @param[in] past trajectoire réalisée jusqu'a la date t
     */
    void asset(PnlMat *path, double t, double T, int nbTimeSteps, PnlRng *rng, const PnlMat *past);

    /**
     * Shift d'une trajectoire du sous-jacent
     *
     * @param[in]  path contient en input la trajectoire
     * du sous-jacent
     * @param[out] shift_path contient la trajectoire path
     * dont la composante d a été shiftée par (1+h)
     * à partir de la date t.
     * @param[in] t date à partir de laquelle on shift
     * @param[in] h pas de différences finies
     * @param[in] d indice du sous-jacent à shifter
     * @param[in] timestep pas de constatation du sous-jacent
     */
    void shiftAsset(PnlMat *shift_path, const PnlMat *path, int d, double h, double t, double timestep);

    /**
     * calcul de Si+u selon le modèle de BS sachant Si et le step u
     *
     * @param[out] double représentant le prix calculé
     * @param[in] trendOrRate qui détermine l'utilisation de la formule
     * avec r_ ou la tendance du modèle
     * @param[in] d l'indice de l'actif permettant d'extraire
     * la ligne convenable de la corrMatrix et la volatilité de
     * cet actif
     */
    double nextPrice(double trendOrRate, double lastPrice, double step, int d, PnlVect *g);

    /**
     * méthode renvoyant une simulation du
     * marché sous la probabilité historique)
     *
     * @param[out] contient la trajectoire du modèle du marché
     * @param[in] nbTimeSteps nombre de pas de constatation
     * @param[in] T date jusqu'à laquelle on simule la trajectoire
     */
    void simul_market(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);
};