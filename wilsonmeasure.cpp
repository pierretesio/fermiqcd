#include "fermiqcd.h"
#include "MyHeaders.h"

int main(int argc, char** argv) {
    mdp.open_wormholes(argc, argv);

    // Paramètres du réseau
    int box[] = {12, 6, 6, 12};   // temps, x, y, z
    mdp_lattice lattice(4, box);

    // Champ de jauge SU(2)
    gauge_field U(lattice, 2);
    coefficients gauge;
    gauge["beta"] = 2.4;

    // Configuration initiale chaude
    set_hot(U);

    // Paramètres du calcul
    int ntherm = 200;     // nombre de sweeps de thermalisation
    int ncfg   = 100;     // nombre de configurations mesurées
    int ndecor = 5;       // sweeps entre deux mesures

    // Tailles de boucles qu'on veut mesurer
    int tvalue[] = {1,2,3,4,5};
    int zvalue[] = {1,2,3,4,5};
    int nt = 5;
    int nz = 5;

    ofstream out("wilson_loops.dat");

    // Thermalisation
    for (int k = 0; k < ntherm; k++) {
        WilsonGaugeAction::heatbath(U, gauge, 1);
        if ((k+1) % 50 == 0) {
            cout << "thermalisation: " << (k+1) << "/" << ntherm << endl;
        }
    }

    // Mesures
    for (int iconf = 1; iconf <= ncfg; iconf++) {

        // Décorrélation entre deux mesures
        WilsonGaugeAction::heatbath(U, gauge, ndecor);

        out << "# configuration " << iconf << endl;

        for (int it = 0; it < nt; it++) {
            for (int iz = 0; iz < nz; iz++) {
                int t = tvalue[it];
                int z = zvalue[iz];

                // boucle rectangulaire dans le plan (0,3)
                double w = my_average_loop(U, 0, t, 3, z);

                out << iconf << " "
                    << t << " "
                    << z << " "
                    << w << endl;
            }
        }

        cout << "mesure " << iconf << "/" << ncfg << endl;
    }

    out.close();
    mdp.close_wormholes();
    return 0;
}