#include <iostream>
#include <fstream>
#include "TropicoFrac.h"

void WriteTeX (std::string nameTeX, TropicoMultiSolve& tms)
{
    std::ifstream inTeXf;
    std::ofstream outTeXf;
    std::string command;
    std::string str;

#ifdef __linux
    command = "mv " + nameTeX + " " + nameTeX + ".c";
#else
    size_t pos = 0;
    std::string nameWithoutPath = nameTeX;

    while ((pos = nameWithoutPath.find('\\')) != std::string::npos)
        nameWithoutPath = nameWithoutPath.substr(pos + 1);

    command = "ren " + nameTeX + " " + nameWithoutPath + ".c";
#endif

    std::system(command.c_str());
    
    try
    {
        inTeXf.open(nameTeX + ".c");
    }
    catch(const std::exception& e)
    {
        throw std::invalid_argument("Incorrect name of TeX.");
    }

    try
    {
        outTeXf.open(nameTeX);
    }
    catch(const std::exception& e)
    {
        throw std::out_of_range("Error opening of file.");
    }

    while (std::getline(inTeXf, str) && str.find("begininsert", 0) == std::string::npos)
        outTeXf << str << "\n";

    outTeXf << str << "\n";

    while (std::getline(inTeXf, str) && str.find("endinsert", 0) == std::string::npos);

    MakeTeX(outTeXf, tms);

    outTeXf << str << "\n";

    while (std::getline(inTeXf, str))
        outTeXf << str << "\n";

    inTeXf.close();
    outTeXf.close();

#ifdef __linux
    command = "rm " + nameTeX + ".c";
#else
    command = "del " + nameTeX + ".c";
#endif

    std::system(command.c_str());
}

int main(int argc, char* argv[])
{
    std::istream* in = & std::cin;
    std::ostream* out = & std::cout;
    std::ifstream inf;
    std::ofstream outf;
#ifdef __linux
    std::string nameTeX = "TeX/Tropico.tex";
#else
    std::string nameTeX = "TeX\\Tropico.tex";
#endif

    if (argc >= 2)
    {
        try
        {
            inf.open(argv[1]);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';

            return 1;
        }
        
        in = &inf;
    }

    if (argc >= 3)
    {
        try
        {
            outf.open(argv[2]);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';

            return 1;
        }
        
        out = &outf;
    }

    if (argc >= 4)
        nameTeX = argv[3];

    if (argc >= 5)
    {
        std::cerr << "Expected: ./name_program [output_TeX_file input_file output_file]\n";

        return 1;
    }

    try
    {
        TropicoMultiSolve tms;

        if (argc == 1)
            *out << "Input matrices:\n";

        *in >> tms;

        tms.Solve();

        *out << tms << "\n";

        WriteTeX(nameTeX, tms);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';

        return 1;
    }

    if (argc >= 2)
        inf.close();
    if (argc >= 3)
        outf.close();

    return 0;
}
