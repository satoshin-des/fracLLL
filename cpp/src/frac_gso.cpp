#include "frac_lll.hpp"

#include <NTL/mat_ZZ.h>

void FracLLL::fracGSO()
{
    long i, j, k, h;
    NTL::ZZ sum;
    NTL::mat_ZZ D, d;

    this->m_G = this->basis * NTL::transpose(this->basis);
    this->B[NUM].SetLength(this->n);
    this->B[DEN].SetLength(this->n);
    this->mu[NUM].SetDims(this->n, this->n);
    this->mu[DEN].SetDims(this->n, this->n);
    D.SetDims(this->n, this->n);

    for (i = 0; i < this->n; ++i)
    {
        this->mu[NUM][i][i] = 1;
        for (j = 0; j < this->n; ++j)
        {
            this->mu[DEN][i][j] = 1;
        }
    }

    for (i = 0; i < this->n; ++i)
    {
        for (j = 0; j < this->n; ++j)
        {
            if ((i == 0) && (j == 0))
            {
                D[i][j] = 1;
            }
            else
            {
                d.SetDims(i, i);
                for (k = 0; k < i; ++k)
                {
                    for (h = 0; h < i; ++h)
                    {
                        if (k < j)
                        {
                            d[k][h] = this->m_G[k][h];
                        }
                        else
                        {
                            d[k][h] = this->m_G[k + 1][h];
                        }
                    }
                }
                D[i][j] = NTL::determinant(d);
            }
        }
    }

    for (i = 0; i < this->n; ++i)
    {
        sum = 0;
        for (j = 0; j <= i; ++j)
        {
            for (k = 0; k <= i; ++k)
            {
                if ((j + k) % 2 == 0)
                {
                    sum += D[i][j] * D[i][k] * this->m_G[j][k];
                }
                else
                {
                    sum -= D[i][j] * D[i][k] * this->m_G[j][k];
                }
            }
        }
        this->B[NUM][i] = sum;
        this->B[DEN][i] = D[i][i] * D[i][i];
        this->m_d = NTL::GCD(this->B[NUM][i], this->B[DEN][i]);
        this->B[NUM][i] /= this->m_d;
        this->B[DEN][i] /= this->m_d;

        for (j = 0; j < i; ++j)
        {
            sum = 0;
            for (k = 0; k <= j; ++k)
            {
                if ((j + k) & 1)
                {
                    sum -= D[j][k] * this->m_G[i][k];
                }
                else
                {
                    sum += D[j][k] * this->m_G[i][k];
                }
            }
            this->mu[NUM][i][j] = D[j][j] * sum;

            sum = 0;
            for (k = 0; k <= j; ++k)
            {
                for (h = 0; h <= j; ++h)
                {
                    if ((k + h) & 1)
                    {
                        sum -= D[j][k] * D[j][h] * this->m_G[k][h];
                    }
                    else
                    {
                        sum += D[j][k] * D[j][h] * this->m_G[k][h];
                    }
                }
            }
            this->mu[DEN][i][j] = sum;

            m_d = NTL::GCD(this->mu[NUM][i][j], this->mu[DEN][i][j]);
            this->mu[NUM][i][j] /= this->m_d;
            this->mu[DEN][i][j] /= this->m_d;
        }
    }
}
