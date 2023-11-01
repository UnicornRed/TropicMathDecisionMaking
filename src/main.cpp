#include <iostream>
#include "TropicoFrac.h"

int main()
{
    try
    {
        TropicoMultiSolve tms;

        std::cin >> tms;

        tms.Solve();

        std::cout << tms << "\n";
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }

    return 0;
}