#include "frac_lll.hpp"

#include <fstream>

void FracLLL::svpChallenge(const long dim, const long seed)
{
    char file_name[200];
    sprintf(file_name, "../../svp_challenge_list/svp_challenge_%ld_%ld.txt", dim, seed);

    std::ifstream fin(file_name);

    this->n = dim;
    this->basis.SetDims(dim, dim);
    for (long i = 0, j; i < dim; ++i)
    {
        for (j = 0; j < dim; ++j)
        {
            fin >> this->basis[i][j];
        }
    }
}
