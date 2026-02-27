#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

ofstream outfile1("calcul_Nd.dat");
ofstream outfile2("sample.dat");
ofstream outfile3("walk_2d.dat");

float myrandom(float a, float b)
{
    float range = pow(2.0, 31.0);
    return a + (b - a) * random() / range;
}

float prob(int d, float *x)
{
    float arg = 0.0f;
    for (int i = 0; i < d; i++) arg = arg + x[i] * x[i];
    return exp(-arg);
}

float func(int d, float *x)
{
    float arg = 0.0f;
    for (int i = 0; i < d; i++) arg = arg + x[i] * x[i];
    return arg;
}

int main()
{
    srandom(12345);

    int j, dim = 2;
    float epsi = 1.0f;

    int i, nmarkov = 100000;
    float x[dim], y[dim];
    float store_func[nmarkov + 1];
    float accept = 0.0f;

    for (j = 0; j < dim; j++) x[j] = 0.0f;

    for (i = 1; i <= nmarkov; i++)
    {
        for (j = 0; j < dim; j++)
            y[j] = x[j] + myrandom(-epsi, epsi);

        if (prob(dim, y) / prob(dim, x) > myrandom(0.0f, 1.0f))
        {
            for (j = 0; j < dim; j++) x[j] = y[j];
            accept = accept + 1.0f;
        }

        store_func[i] = func(dim, x);

        if (dim == 1) outfile2 << i << " " << x[0] << endl;
        if (dim == 2) outfile3 << x[0] << " " << x[1] << endl;
    }

    float sum = 0.0f, mean = 0.0f, sum2 = 0.0f, mean2 = 0.0f, sigma = 0.0f;
    float exact = dim * 0.5f;

    for (i = 1; i <= nmarkov; i++)
    {
        sum = sum + store_func[i];
        mean = sum / i;

        sum2 = sum2 + store_func[i] * store_func[i];
        mean2 = sum2 / i;

        sigma = sqrt((mean2 - mean * mean) / i);

        if (i % 1000 == 0)
            outfile1 << i << " " << mean << " " << sigma << endl;
    }

    cout << "calcul en dimension: " << dim << endl;
    cout << "accept rate: " << accept / nmarkov << endl;
    cout << "mean: " << mean << " std error: " << sigma << " exact: " << exact << endl;

    outfile1.close();
    outfile2.close();
    outfile3.close();

    return 0;
}