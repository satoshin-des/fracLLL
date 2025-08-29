#include <iostream>

#include "frac_lll.hpp"

int main(int argc, const char** argv)
{
    FracLLL L;

    L.svpChallenge(atol(argv[1]), 0);
    L.fracDeepLLL(99, 100);

    std::cout << L.basis << std::endl;

    return 0;
}
