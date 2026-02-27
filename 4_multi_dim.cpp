#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

float myrandom(float a, float b)
{
    float range = pow(2.0, 31.0);
    return a + (b - a) * random() / range;
}

float prob(int d, float *x)
{
    float arg = 0.0f;
    for (int i = 0; i < d; i++) arg += x[i] * x[i];
    return exp(-arg);
}

float func(int d, float *x)
{
    float arg = 0.0f;
    for (int i = 0; i < d; i++) arg += x[i] * x[i];
    return arg;
}

int main()
{
    srandom(12345);

    float epsi = 2.0f;
    int nmarkov = 100000;
    int ncorel = 100;

    for (int dim = 1; dim <= 5; dim++)
    {
        string f_calcul = "calcul_Nd_D" + to_string(dim) + ".dat";
        string f_sample = "sample_D" + to_string(dim) + ".dat";
        string f_walk2d  = "walk_2d_D" + to_string(dim) + ".dat";
        string f_autoc   = "autocorel_Nd_D" + to_string(dim) + ".dat";

        ofstream outfile1(f_calcul.c_str());
        ofstream outfile2(f_sample.c_str());
        ofstream outfile3(f_walk2d.c_str());
        ofstream outfile4(f_autoc.c_str());

        float *x = new float[dim];
        float *y = new float[dim];
        float *store_func = new float[nmarkov + 1];

        for (int j = 0; j < dim; j++) x[j] = 0.0f;

        float accept = 0.0f;

        for (int i = 1; i <= nmarkov; i++)
        {
            for (int j = 0; j < dim; j++)
                y[j] = x[j] + myrandom(-epsi, epsi);

            if (prob(dim, y) / prob(dim, x) > myrandom(0.0f, 1.0f))
            {
                for (int j = 0; j < dim; j++) x[j] = y[j];
                accept += 1.0f;
            }

            store_func[i] = func(dim, x);

            if (dim == 1) outfile2 << i << " " << x[0] << endl;
            if (dim == 2) outfile3 << x[0] << " " << x[1] << endl;
        }

        float sum = 0.0f, mean = 0.0f, sum2 = 0.0f, mean2 = 0.0f, sigma = 0.0f;
        float exact = dim * 0.5f;

        for (int i = 1; i <= nmarkov; i++)
        {
            sum += store_func[i];
            mean = sum / i;

            sum2 += store_func[i] * store_func[i];
            mean2 = sum2 / i;

            sigma = sqrt((mean2 - mean * mean) / i);

            if (i % 1000 == 0)
                outfile1 << i << " " << mean << " " << sigma << endl;
        }

        float sumshift, chi, chi0 = 0.0f, normedchi, rough_time = 0.0f;

        for (int icorel = 0; icorel <= ncorel; icorel++)
        {
            sumshift = 0.0f;
            for (int i = 1; i <= nmarkov - icorel; i++)
                sumshift += store_func[i] * store_func[i + icorel];

            chi = sumshift / (nmarkov - icorel) - mean * mean;
            if (icorel == 0) chi0 = chi;
            normedchi = chi / chi0;

            if (icorel == 1 && normedchi > 0.0f)
                rough_time = -1.0f / log(normedchi);

            outfile4 << icorel << " " << normedchi << endl;
        }

        cout << "dimension: " << dim << endl;
        cout << "accept rate: " << accept / nmarkov << endl;
        cout << "mean: " << mean << " std error: " << sigma << " exact: " << exact << endl;
        cout << "rough estimate of exponential time: " << rough_time << endl;
        cout << endl;

        outfile1.close();
        outfile2.close();
        outfile3.close();
        outfile4.close();

        delete[] x;
        delete[] y;
        delete[] store_func;
    }

    return 0;
}