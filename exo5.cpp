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

    int dim = 1;
    float epsi = 2.0f;

    int nsample = 100000;
    int nskip = 0;

    int ncorel = 100;

    int minBlock = 10;
    int maxBlock = 5000;
    int stepBlock = 10;

    float *x = new float[dim];
    float *y = new float[dim];
    float *store_func = new float[nsample + 1];

    for (int j = 0; j < dim; j++) x[j] = 0.0f;

    float accept = 0.0f;
    int ntry = 0;

    for (int t = 1; t <= nsample; t++)
    {
        for (int s = 0; s <= nskip; s++)
        {
            for (int j = 0; j < dim; j++)
                y[j] = x[j] + myrandom(-epsi, epsi);

            ntry++;

            if (prob(dim, y) / prob(dim, x) > myrandom(0.0f, 1.0f))
            {
                for (int j = 0; j < dim; j++) x[j] = y[j];
                accept += 1.0f;
            }
        }

        store_func[t] = func(dim, x);

        if (dim == 1) outfile2 << t << " " << x[0] << endl;
        if (dim == 2) outfile3 << x[0] << " " << x[1] << endl;
    }

    float sum = 0.0f, mean = 0.0f, sum2 = 0.0f, mean2 = 0.0f, sigma = 0.0f;
    float exact = dim * 0.5f;

    for (int i = 1; i <= nsample; i++)
    {
        sum += store_func[i];
        mean = sum / i;

        sum2 += store_func[i] * store_func[i];
        mean2 = sum2 / i;

        sigma = sqrt((mean2 - mean * mean) / i);

        if (i % 1000 == 0)
            outfile1 << i << " " << mean << " " << sigma << endl;
    }

    cout << "dimension: " << dim << endl;
    cout << "epsi: " << epsi << endl;
    cout << "nsample: " << nsample << " nskip: " << nskip << endl;
    cout << "accept rate (per proposal): " << accept / ntry << endl;
    cout << "mean: " << mean << " std error: " << sigma << " exact: " << exact << endl;

    ofstream outA("autocorel_Nd.dat");
    float sumshift,chi,chi0,normedchi,rough_time=0.0f;
    for(int icorel=0;icorel<=ncorel;icorel++)
    {
        sumshift=0.0f;
        for(int i=1;i<=nsample-icorel;i++)
            sumshift += store_func[i]*store_func[i+icorel];

        chi = sumshift/(nsample-icorel) - mean*mean;
        if(icorel==0) chi0=chi;
        normedchi = chi/chi0;

        if(icorel==1 && normedchi>0.0f)
            rough_time = -1.0f/log(normedchi);

        outA << icorel << " " << normedchi << endl;
    }
    cout << "rough estimate of exponential time: " << rough_time << endl;
    outA.close();

    ofstream outJK("jackknife.dat");

    float totalSum = 0.0f;
    for (int i = 1; i <= nsample; i++) totalSum += store_func[i];

    for (int B = minBlock; B <= maxBlock; B += stepBlock)
    {
        int nb = nsample / B;
        if (nb < 2) break;

        float *blockSum = new float[nb];
        for (int k = 0; k < nb; k++) blockSum[k] = 0.0f;

        for (int k = 0; k < nb; k++)
        {
            int start = k * B + 1;
            int end = (k + 1) * B;
            float s = 0.0f;
            for (int i = start; i <= end; i++) s += store_func[i];
            blockSum[k] = s;
        }

        float meanJK = 0.0f;
        float *mleave = new float[nb];

        for (int k = 0; k < nb; k++)
        {
            float m = (totalSum - blockSum[k]) / (nsample - B);
            mleave[k] = m;
            meanJK += m;
        }
        meanJK /= nb;

        float v = 0.0f;
        for (int k = 0; k < nb; k++)
        {
            float d = mleave[k] - meanJK;
            v += d * d;
        }
        v *= (nb - 1.0f) / nb;

        float errJK = sqrt(v);

        outJK << B << " " << errJK << " " << nb << endl;

        delete[] blockSum;
        delete[] mleave;
    }

    outJK.close();

    outfile1.close();
    outfile2.close();
    outfile3.close();

    delete[] x;
    delete[] y;
    delete[] store_func;

    return 0;
}