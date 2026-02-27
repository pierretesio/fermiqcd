// PROGRAM1: Use and test of random generator + calcul des integrales (Monte-Carlo)

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>   // random(), srandom()

using namespace std;

ofstream outfile("random.dat");   // open stream for output

// random number in [a,b]
float myrandom(float a, float b)
{
    float range = pow(2.0, 31.0);
    return a + (b - a) * random() / range;
}

int main()
{
    // IMPORTANT: initialiser le generateur (sinon reproductible "par defaut")
    srandom(12345);

    // ---- Partie 1: test generateur + fichier pour histogramme ----
    float x, sum = 0.0f, mean;
    int i, hitnb = 100000;        // make hitnb calls

    cout << "create " << hitnb << " rand. numb. in [0,1]\n"
         << "and compute mean" << endl;

    for (i = 1; i <= hitnb; i++)
    {
        x = myrandom(0.0f, 1.0f);
        outfile << i << " " << x << endl;   // writ[e] file for histogram

        sum = sum + x;
        mean = sum / i;

        if (i % (hitnb / 100) == 0)
        {
            cout << "i= " << i << " mean= " << mean << endl;
        }
    }

    outfile.close();

    // ---- Partie 2: calcul des integrales I0 et I2 par Monte-Carlo ----
    // I_n = int_0^1 exp(-x^2) x^n dx, n=0,2
    int N = 1000000; // tu peux changer
    float s0 = 0.0f, s2 = 0.0f;

    for (i = 1; i <= N; i++)
    {
        x = myrandom(0.0f, 1.0f);
        float fx = exp(-x * x);
        s0 += fx;           // n=0
        s2 += fx * x * x;   // n=2
    }

    float I0 = s0 / N;
    float I2 = s2 / N;

    cout << "\nMC integrals with N=" << N << endl;
    cout << "I0 ~ " << I0 << " (exact 0.746824)" << endl;
    cout << "I2 ~ " << I2 << " (exact 0.189472)" << endl;

    return 0;
}