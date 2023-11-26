#include <exception>
#include <sstream>
#include <map>
#include <iomanip>
#include <numeric>
#include "TropicoFrac.h"

void GcdFrac(int& n, int& m)
{
    int minNM = std::min(std::abs(n), m), maxDel = 1;

    for (int i = 2; i <= minNM; ++i)
        if (n % i == 0 && m % i == 0)
            maxDel = i;

    n /= maxDel;
    m /= maxDel;

    if (n == 0)
        m = 1;
}

void DegreeFrac::Reduction()
{
    if (m < 0)
    {
        m *= -1;
        n *= -1;
    }

    GcdFrac(n, m);
}

DegreeFrac& DegreeFrac::operator+=(const DegreeFrac& df)
{
    n = n * df.m + df.n * m;
    m *= df.m;

    Reduction();

    return *this;
}

DegreeFrac DegreeFrac::operator+(const DegreeFrac& df) const
{
    DegreeFrac res(*this);

    return res += df;
}

DegreeFrac& DegreeFrac::operator*=(const DegreeFrac& df)
{
    n *= df.n;
    m *= df.m;

    Reduction();

    return *this;
}

DegreeFrac DegreeFrac::operator*(const DegreeFrac& df) const
{
    DegreeFrac res(*this);

    return res *= df;
}

DegreeFrac DegreeFrac::operator-() const
{
    DegreeFrac res(*this);

    res.n *= -1;

    return res;
}

bool DegreeFrac::operator<=(const DegreeFrac& df) const
{
    return n * df.m <= df.n * m;
}

std::ostream& operator<<(std::ostream& out, const DegreeFrac& df)
{
    out << "(" << df.n;

    if (df.m != 1)
        out << "/" << df.m;
        
    out << ")";

    return out;
}

std::istream& operator>>(std::istream& in, DegreeFrac& df)
{
    size_t pos;
    std::string str;

    in >> str;

    try
    {
        pos = str.find('/');
        df.n = std::stoi(str.substr(0, pos));
        str = str.substr(pos + 1);

        if (pos != std::string::npos)
            df.m = std::stoi(str);
    }
    catch(const std::exception& e)
    {
        throw std::invalid_argument("Invalid format. Try: int[/int]. For example: 1/2.");
    }

    df.Reduction();

    return in;
}

PrimeDegreeFrac& PrimeDegreeFrac::operator*=(const PrimeDegreeFrac& pdf)
{
    df += pdf.df;

    return *this;
}

PrimeDegreeFrac PrimeDegreeFrac::operator*(const PrimeDegreeFrac& pdf) const
{
    PrimeDegreeFrac res(*this);

    return res *= pdf;
}

PrimeDegreeFrac& PrimeDegreeFrac::operator^=(const DegreeFrac& _df)
{
    df *= _df;

    return *this;
}

PrimeDegreeFrac PrimeDegreeFrac::operator^(const DegreeFrac& df) const
{
    PrimeDegreeFrac res(*this);

    return res^df;
}

bool PrimeDegreeFrac::operator<=(const PrimeDegreeFrac& pdf) const
{
    return n <= pdf.n;
}

TropicoFrac::TropicoFrac(int _n, int _m)
{
    if (_m == 0)
        throw std::invalid_argument("The denominator cannot be 0.");

    if (_n < 0 || _m < 0)
        throw std::invalid_argument("The number cannot be less 0.");

    GcdFrac(_n, _m);
    Convert(_n, 1);
    Convert(_m, -1);
}

bool TropicoFrac::ContainPrime(int _n) const
{
    for (const PrimeDegreeFrac& p : n)
    {
        if (p.GetPrime() == _n)
            return true;
    }

    return false;
}

DegreeFrac TropicoFrac::GetDegreePrime(int _n) const
{
    DegreeFrac df;

    for (const PrimeDegreeFrac& p : n)
    {
        if (p.GetPrime() == _n)
            return p.GetDegreeFrac();
    }

    throw std::out_of_range("It not have this prime number.");
}

void TropicoFrac::AddDegree(int prime, DegreeFrac df)
{
    if (prime == 1)
        return;

    if (ContainPrime(prime))
    {
        DegreeFrac EraseDf = GetDegreePrime(prime);
        df += EraseDf;
        n.erase(PrimeDegreeFrac(prime, EraseDf));
    }

    if (df != DegreeFrac(0, 1))
        n.insert(PrimeDegreeFrac(prime, df));
}

void TropicoFrac::Convert(int _n, const DegreeFrac& delta)
{
    if (_n == 0)
        AddDegree(0, 1);

    for (int i = 2; i <= _n; ++i)
        while (_n % i == 0)
        {
            AddDegree(i, delta);
            _n /= i;
        }
}

void TropicoFrac::PositiveDegree(TropicoFrac& tf1, TropicoFrac& tf2) const
{
    DegreeFrac df;

    for (const PrimeDegreeFrac& p : n)
    {
        df = p.GetDegreeFrac();

        if (df < DegreeFrac(0, 1))
        {
            tf1.AddDegree(p.GetPrime(), -df);
            tf2.AddDegree(p.GetPrime(), -df);
        }
    }
}

bigint TropicoFrac::ValueIntPositive() const
{
    bigint res(1);

    if (ContainPrime(0))
        return bigint(0);

    for (const PrimeDegreeFrac& p : n)
        for (int i{}; i < p.GetDegreeFrac().GetN(); ++i)
            res *= p.GetPrime();

    return res;
}

void SqrtProd(std::ostream& out, const std::pair<int, int> t, int sign = 1)
{
    if (t.first != sign)
    {
        if (std::abs(t.first) == 2)
            out << "\\sqrt";
        else
            out << "\\sqrt[" << std::abs(t.first) << "]";
    }
    
    out << "{" << t.second << "}";
}

double TropicoFrac::ToDouble() const
{
    double res = 1;

    for (const PrimeDegreeFrac& p : n)
        res *= p.ToDouble();

    return res;
}

std::ostream& MakeTeX(std::ostream& out, const TropicoFrac& tf)
{
    bool isFrac = false, positFrac = false;
    int signN, signM;
    std::map<int, int> TeXPDF;

    for (const PrimeDegreeFrac& p : tf.n)
    {
        if (p.GetDegreeFrac() < 0)
            isFrac = true;
        
        signN = p.GetDegreeFrac().GetN() >= 0 ? 1 : -1;
        signM = signN * p.GetDegreeFrac().GetM();

        if (!TeXPDF.contains(signM))
            TeXPDF[signM] = 1;

        for (int i{}; i < signN * p.GetDegreeFrac().GetN(); ++i)
            TeXPDF[signM] *= p.GetPrime();
    }

    out << "{";

    for (const auto& t : TeXPDF)
        if (t.first > 0)
        {
            positFrac = true;
            SqrtProd(out, t);
        }

    if (!positFrac)
        out << 1;

    out << "}";

    if (isFrac)
        out << " / ";

    out << "{";

    for (auto t = TeXPDF.crbegin(); t != TeXPDF.crend(); ++t)
        if (t->first < 0)
            SqrtProd(out, *t, -1);

    out << "}";

    return out;
}

bool TropicoFrac::operator<=(const TropicoFrac& tf) const
{
    DegreeFrac df;
    TropicoFrac tf1(*this), tf2(tf), tf3;

    this->PositiveDegree(tf1, tf2);
    tf.PositiveDegree(tf1, tf2);

    int m = 1;

    for (const PrimeDegreeFrac& p : n)
        m = std::lcm(m, p.GetDegreeFrac().GetM());

    for (const PrimeDegreeFrac& p : tf.n)
        m = std::lcm(m, p.GetDegreeFrac().GetM());

    tf1 ^= DegreeFrac(m);
    tf2 ^= DegreeFrac(m);
    tf3 = tf1;

    for (const PrimeDegreeFrac& p : tf3.n)
    {
        if (tf2.ContainPrime(p.GetPrime()))
        {
            df = tf1.GetDegreePrime(p.GetPrime()) <= tf2.GetDegreePrime(p.GetPrime()) ?
                 tf1.GetDegreePrime(p.GetPrime()) : tf2.GetDegreePrime(p.GetPrime());
            tf1 *= PrimeDegreeFrac(p.GetPrime(), -df);
            tf2 *= PrimeDegreeFrac(p.GetPrime(), -df);
        }
    }
    bigint val1(tf1.ValueIntPositive()), val2(tf2.ValueIntPositive());
    bool res = val1 <= val2;
    if (res != (this->ToDouble() <= tf.ToDouble()))
    {
        std::cout << "ERROR!!!\n";
        std::cout << m << "\n";
        std::cout << tf1 << " " << tf2 << "\n";
        std::cout << tf1.ValueIntPositive() << " " << tf2.ValueIntPositive() << "\n";
        std::cout << tf1.ToDouble() << " " << tf2.ToDouble() << "\n";
        std::cout << *this << " " << tf << "\n";
        std::cout << this->ToDouble() << " " << tf.ToDouble() << "\n\n";
    }

    return res;
}

TropicoFrac& TropicoFrac::operator+=(const TropicoFrac& tf)
{
    *this = *this >= tf ? *this : tf;

    return *this;
}

TropicoFrac TropicoFrac::operator+(const TropicoFrac& tf) const
{
    TropicoFrac res(*this);

    return res += tf;
}

TropicoFrac& TropicoFrac::operator*=(const TropicoFrac& tf)
{
    TropicoFrac res(*this);

    if (tf.ContainPrime(0))
    {
        n.clear();
        n.insert(PrimeDegreeFrac(0, DegreeFrac()));

        return *this;
    }

    if (ContainPrime(0))
        return *this;

    for (const PrimeDegreeFrac& p : tf.n)
        res.AddDegree(p.GetPrime(), p.GetDegreeFrac());

    this->n = res.n;

    return *this;
}

TropicoFrac TropicoFrac::operator*(const TropicoFrac& tf) const
{
    TropicoFrac res(*this);

    return res *= tf;
}

TropicoFrac& TropicoFrac::operator^=(const DegreeFrac& df)
{
    std::set<PrimeDegreeFrac> res;

    for (const PrimeDegreeFrac& p : n)
        res.insert(PrimeDegreeFrac(p.GetPrime(), p.GetDegreeFrac() * df));
    
    n = res;

    return *this;
}

TropicoFrac TropicoFrac::operator^(const DegreeFrac& df) const
{
    TropicoFrac res(*this);

    return res ^= df;
}

std::ostream& operator<<(std::ostream& out, const TropicoFrac& tf)
{
    bool first = true;

    if (!tf.n.size())
        out << 1;

    for (const PrimeDegreeFrac& p : tf.n)
    {
        if (first)
            first = false;
        else
            out << "*";

        out << p.GetPrime();

        if (p.GetDegreeFrac() != DegreeFrac(1, 1))
            out << "^" << p.GetDegreeFrac();
    }

    return out;
}

std::istream& operator>>(std::istream& in, TropicoFrac& tf)
{
    size_t pos;
    int n, m;
    DegreeFrac df;
    std::stringstream ss;
    std::string str, token;

    tf.n.clear();

    in >> str;

    try
    {
        if ((pos = str.find('/')) == std::string::npos && str.find('^') == std::string::npos)
        {
            n = std::stoi(str.substr(0));

            tf = TropicoFrac(n, 1);
        }
        else if (str.find('^') == std::string::npos)
        {
            n = std::stoi(str.substr(0, pos));
            m = std::stoi(str.substr(pos + 1));

            tf = TropicoFrac(n, m);
        }
        else
        {
        do
            {
                df = DegreeFrac();

                pos = str.find('*');
                token = str.substr(0, pos);

                if (pos != std::string::npos)
                    str = str.substr(pos + 1);
                else
                    str = "";

                pos = token.find('^');
                n = std::stoi(token.substr(0, pos));
                token = token.substr(pos + 1);
                ss.flush();

                if (pos != std::string::npos)
                {
                    ss.clear();
                    ss << token.substr(0);
                    ss >> df;
                }

                tf.AddDegree(n, df);
            } while (!str.empty());
        }
    }
    catch(const std::exception& e)
    {
        throw std::invalid_argument("Invalid format. Try: int[^int[/int]][*...*][int[^int[/int]]] or int/int. For example: 2^3/4*5^6/7.");
    }

    return in;
}

TropicoMatrix::TropicoMatrix(size_t height, size_t width, TropicoFrac elem)
{
    matr.resize(height, std::vector<TropicoFrac>(width, elem));
}

size_t TropicoMatrix::GetWidth() const
{
    if (GetHeight() == 0)
        return 0;
    else
        return matr[0].size();
}

TropicoMatrix TropicoMatrix::GetI(size_t size)
{
    TropicoMatrix res(size, size);

    for (size_t i{}; i < size; ++i)
        res[i][i] = 1;

    return res;
}

void TropicoMatrix::Resize(size_t height, size_t width)
{
    matr.resize(height);

    std::for_each(matr.begin(), matr.end(), [width](std::vector<TropicoFrac>& v)
    {
        v.resize(width, TropicoFrac());
    });
}

TropicoMatrix TropicoMatrix::SubMatrix(size_t x, size_t y, size_t height, size_t width) const
{
    if (x + height > GetHeight() || y + width > GetWidth())
        throw std::invalid_argument("SubMatrix must be less than Matrix.");

    TropicoMatrix sub(height, width);

    for (size_t i{}; i < height; ++i)
        for (size_t j{}; j < width; ++j)
            sub[i][j] = matr[i + x][j + y];

    return sub;
}

TropicoMatrix TropicoMatrix::MultiConjTrans() const
{
    TropicoMatrix MCT(GetWidth(), GetHeight());

    for (size_t i{}; i < MCT.GetHeight(); ++i)
        for (size_t j{}; j < MCT.GetWidth(); ++j)
            if (matr[j][i] != 0)
                MCT[i][j] = matr[j][i] ^ (-1);
            else
                MCT[i][j] = 0;

    return MCT;
}

TropicoMatrix& TropicoMatrix::Standardization()
{
    if (GetHeight() < 1 || GetWidth() < 1)
        throw std::invalid_argument("Matrix is empty.");

    for (size_t j{}; j < GetWidth(); ++j)
    {
        TropicoFrac sum;
        TropicoMatrix _tm = SubMatrix(0, j, GetHeight(), 1);

        for (size_t i{}; i < GetHeight(); ++i)
            sum += matr[i][j];
        
        _tm *= sum ^ (-1);

        for (size_t i{}; i < GetHeight(); ++i)
            matr[i][j] = _tm[i][0];
    }

    return *this;
}

std::vector<std::vector<double>> TropicoMatrix::ToDouble() const
{
    std::vector<std::vector<double>> res(GetHeight(), std::vector<double>(GetWidth()));

    for (size_t i{}; i < GetHeight(); ++i)
        for(size_t j{}; j < GetWidth(); ++j)
            res[i][j] = matr[i][j].ToDouble();

    return res;
}

std::ostream& ToDoubleTeXOut(std::ostream& out, const TropicoMatrix& tm)
{
    std::vector<std::vector<double>> outMatrix = tm.ToDouble();

    out << "\t\t\\begin{pmatrix}\n";

    for (size_t i{}; i < outMatrix.size(); ++i)
    {
        out << "\t\t\t";

        for(size_t j{}; j < outMatrix[0].size(); ++j)
            out << std::setprecision(4) << std::setw(8) << std::left << double(outMatrix[i][j]) << (j == tm.GetWidth() - 1 ? " " : " & ");

        out << "\\\\\n";
    }

    out << "\t\t\\end{pmatrix}";

    return out;
}

std::ostream& ToDoubleOut(std::ostream& out, const TropicoMatrix& tm)
{
    std::vector<std::vector<double>> outMatrix = tm.ToDouble();

    for (size_t i{}; i < outMatrix.size(); ++i)
    {
        for(size_t j{}; j < outMatrix[0].size(); ++j)
            out << std::setprecision(4) << std::setw(8) << std::left << double(outMatrix[i][j]) << " ";

        out << "\n";
    }

    return out;
}

std::ostream& MakeTeX(std::ostream& out, const TropicoMatrix& tm)
{
    out << "\t\t\\begin{pmatrix}\n";

    for (size_t i{}; i < tm.GetHeight(); ++i)
    {
        out << "\t\t\t";

        for(size_t j{}; j < tm.GetWidth(); ++j)
            MakeTeX(out,tm.matr[i][j]) << (j == tm.GetWidth() - 1 ? " " : " & ");

        out << "\\\\\n";
    }

    out << "\t\t\\end{pmatrix}";

    return out;
}

bool CorrelVector(const TropicoMatrix& tm1, const TropicoMatrix& tm2)
{
    if (tm1.GetWidth() != 1 || tm2.GetWidth() != 1)
        throw std::invalid_argument("It is not vector.");

    if (tm1.GetHeight() < 1 || tm2.GetHeight() < 1)
        throw std::invalid_argument("Vector is empty.");

    if (tm1.GetHeight() != tm2.GetHeight())
        throw std::invalid_argument("Size of vectors are not equal.");

    TropicoFrac coefProp(tm1[0][0] * (tm2[0][0] ^ (-1)));

    for (size_t i{}; i < tm1.GetHeight(); ++i)
        if (tm1[i][0] * (tm2[i][0] ^ (-1)) != coefProp)
            return false;

    return true;
}

TropicoMatrix& TropicoMatrix::operator+=(const TropicoMatrix& tm)
{
    if (GetHeight() != tm.GetHeight() || GetWidth() != tm.GetWidth())
        throw std::invalid_argument("Size of matrixes not equal.");

    for (size_t i{}; i < tm.GetHeight(); ++i)
        for(size_t j{}; j < tm.GetWidth(); ++j)
            matr[i][j] += tm[i][j];

    return *this;
}

TropicoMatrix TropicoMatrix::operator+(const TropicoMatrix& tm) const
{
    TropicoMatrix res(*this);

    return res += tm;
}

TropicoMatrix& TropicoMatrix::operator*=(const TropicoMatrix& tm)
{
    *this = *this * tm;

    return *this;
}

TropicoMatrix TropicoMatrix::operator*(const TropicoMatrix& tm) const
{
    if (GetWidth() != tm.GetHeight())
        throw std::invalid_argument("Width of first matrix and height of second matrix not equal.");

    TropicoMatrix res(GetHeight(), tm.GetWidth());

    for (size_t i{}; i < res.GetHeight(); ++i)
        for(size_t j{}; j < res.GetWidth(); ++j)
            for (size_t l{}; l < GetWidth(); ++l)
                res[i][j] += matr[i][l] * tm[l][j];
        
    return res;
}

TropicoMatrix& TropicoMatrix::operator*=(const TropicoFrac& tf)
{
    for (size_t i{}; i < GetHeight(); ++i)
        for(size_t j{}; j < GetWidth(); ++j)
            matr[i][j] *= tf;

    return *this;
}

TropicoMatrix TropicoMatrix::operator*(const TropicoFrac& tf) const
{
    TropicoMatrix res(*this);

    return res *= tf;
}

TropicoFrac TropicoMatrix::Trace() const
{
    if (GetWidth() != GetHeight())
        throw std::invalid_argument("Width of matrix and height of matrix not equal.");

    TropicoFrac res;

    for (size_t i{}; i < GetHeight(); ++i)
        res += matr[i][i];

    return res;
}

TropicoFrac TropicoMatrix::SpectralRadius(std::vector<TropicoMatrix>* degreeMatrix) const
{
    if (GetWidth() != GetHeight())
        throw std::invalid_argument("Width of matrix and height of matrix not equal.");

    TropicoFrac res;
    TropicoMatrix C(*this);

    for (size_t i{}; i < GetHeight(); ++i, C *= *this)
    {
        if (degreeMatrix)
            degreeMatrix->push_back(C);
        
        res += C.Trace() ^ DegreeFrac(1, i + 1);
    }

    return res;
}

TropicoMatrix TropicoMatrix::Kleene() const
{
    if (GetWidth() != GetHeight())
        throw std::out_of_range("Width of matrix and height of matrix not equal.");

    TropicoMatrix res = TropicoMatrix::GetI(GetHeight());
    TropicoMatrix C(*this);

    for (size_t i{1}; i < GetHeight(); ++i, C *= *this)
        res += C;

    return res;
}

TropicoMatrix TropicoMatrix::RemoveCorrel() const
{
    bool isCorrel = false;
    TropicoMatrix _tm(SubMatrix(0, 0, GetHeight(), 1));

    for (size_t i{1}; i < GetWidth(); ++i)
    {
        isCorrel = false;

        for (size_t j{}; j < _tm.GetWidth(); ++j)
            if (CorrelVector(SubMatrix(0, i, GetHeight(), 1), _tm.SubMatrix(0, j, _tm.GetHeight(), 1)))
                isCorrel = true;

        if (!isCorrel)
        {
            _tm.Resize(_tm.GetHeight(), _tm.GetWidth() + 1);

            for (size_t j{}; j < _tm.GetHeight(); ++j)
                _tm[j][_tm.GetWidth() - 1] = matr[j][i];
        }
    }

    return _tm;
}

std::ostream& operator<<(std::ostream& out, const TropicoMatrix& tm)
{
    size_t maxLength = 0;
    std::string str;
    std::stringstream ss;
    std::vector<std::string> strMatrix;

    for (size_t i{}; i < tm.GetHeight(); ++i)
    {
        for(size_t j{}; j < tm.GetWidth(); ++j)
        {
            ss.clear();
            ss << tm[i][j];
            std::getline(ss, str);

            maxLength = std::max(maxLength, str.length());
            
            strMatrix.push_back(str);
        }
    }

    for (size_t i{}; i < tm.GetHeight(); ++i)
    {
        for(size_t j{}; j < tm.GetWidth(); ++j)
        {
            strMatrix[j + i * tm.GetWidth()].resize(maxLength, ' ');

            out << strMatrix[j + i * tm.GetWidth()] << " ";
        }
        
        out << std::endl;
    }

    return out;
}

std::istream& operator>>(std::istream& in, TropicoMatrix& tm)
{
    size_t height, width;

    in >> height >> width;

    tm.Resize(height, width);

    for (size_t i{}; i < tm.GetHeight(); ++i)
        for(size_t j{}; j < tm.GetWidth(); ++j)
            in >> tm[i][j];

    return in;
}

TropicoSolve::TropicoSolve(TropicoMatrix _tm) : tm(_tm)
{
    if (tm.GetWidth() != tm.GetHeight())
        throw std::invalid_argument("Width of matrix and height of matrix not equal.");
}

std::vector<TropicoMatrix>& TropicoSolve::WorstSolve()
{
    delta = (TropicoMatrix(1, tm.GetWidth(), 1) * kleene * TropicoMatrix(tm.GetHeight(), 1, 1))[0][0];
    worstMatr = ((delta ^ (-1)) * TropicoMatrix(tm.GetHeight(), tm.GetWidth(), 1) + (spectral ^ (-1)) * tm).Kleene();

    TropicoMatrix worstMatrWithoutCor = worstMatr.RemoveCorrel();
    
    for (size_t i{}; i < worstMatrWithoutCor.GetWidth(); ++i)
        worst.push_back((worstMatrWithoutCor.SubMatrix(0, i, worstMatrWithoutCor.GetHeight(), 1) * TropicoFrac(1)).Standardization());

    return worst;
}

TropicoMatrix TropicoSolve::MakeP() const
{
    TropicoMatrix _P;

    _P = kleene.RemoveCorrel();

    return _P;
}

std::vector<TropicoMatrix> TropicoSolve::MakePlk() const
{
    std::vector<size_t> k, l;
    TropicoFrac maxPj, value;
    std::vector<TropicoMatrix> _Plk;

    for (size_t i{}; i < P.GetWidth(); ++i)
    {
        if ((value = (TropicoMatrix(1, tm.GetWidth(), 1) * (P.SubMatrix(0, i, P.GetHeight(), 1) * P.SubMatrix(0, i, P.GetHeight(), 1).MultiConjTrans()) *
                      TropicoMatrix(tm.GetHeight(), 1, 1))[0][0]) >= maxPj)
        {
            if (value == maxPj)
                k.push_back(i);
            else
            {
                k.clear();
                k.push_back(i);
            }

            maxPj = value;    
        }
    }

    maxPj = 0;

    for (size_t i{}; i < P.GetHeight(); ++i)
    {
        for (const size_t& _k : k)
            if ((value = P[i][_k] ^ (-1)) > maxPj)
            {
                if (value == maxPj)
                    l.push_back(i);
                else
                {
                    l.clear();
                    l.push_back(i);
                }

                maxPj = value;
            }
    }

    for (const size_t& _k : k)
    {
        for (const size_t& _l : l)
        {
            TropicoMatrix _P(P.GetHeight(), P.GetWidth());
            _P[_l][_k] = P[_l][_k];
            _Plk.push_back(_P);
        }
    }

    return _Plk;
}

std::vector<std::vector<TropicoMatrix>>& TropicoSolve::BestSolve()
{
    P = MakeP();
    Plk = MakePlk();

    for (const TropicoMatrix& _tm : Plk)
    {
        TropicoMatrix bestMatrWithoutCor = P * (TropicoMatrix::GetI(P.GetWidth()) + _tm.MultiConjTrans() * P);
        bestMatr.push_back(bestMatrWithoutCor);
        bestMatrWithoutCor = bestMatrWithoutCor.RemoveCorrel();
        best.push_back(std::vector<TropicoMatrix>());

        for (size_t i{}; i < bestMatrWithoutCor.GetWidth(); ++i)
            best.back().push_back((bestMatrWithoutCor.SubMatrix(0, i, bestMatrWithoutCor.GetHeight(), 1) * TropicoFrac(1)).Standardization());
    }

    return best;
}

void TropicoSolve::Solve()
{
    spectral = tm.SpectralRadius(&degreeTM);
    kleene = ((spectral ^ (-1)) * tm).Kleene();
    kleeneWithoutCorrel = kleene.RemoveCorrel();
    kleeneWithoutCorrel.Standardization();

    WorstSolve();
    BestSolve();
}

std::ostream& MakeTeXMatrix(std::ostream& out, std::string str, const TropicoMatrix& tm)
{
    out << "\t\\[\n";
    out << "\t\t" << str << " = \n";
    MakeTeX(out, tm) << "\n\t\t\\approx\n";
    ToDoubleTeXOut(out, tm) << "\n\t\\]\n";

    return out;
}

std::ostream& MakeTeXFrac(std::ostream& out, std::string str, const TropicoFrac& tf)
{
    out << "\t\\[\n";
    out << "\t\t" << str << " = ";
    MakeTeX(out, tf) << "\\approx " << tf.ToDouble() << "\n\t\\]\n";

    return out;
}

std::ostream& MakeTeX(std::ostream& out, const TropicoSolve& ts)
{
    MakeTeXMatrix(out, "\\bold{" + ts.GetName() + "}", ts.GetMatrix());

    for (size_t i{1}; i < ts.GetDegreeMatrix().size(); ++i)
        MakeTeXMatrix(out, "\\left(\\bold{" + ts.GetName() + "}\\right)^{" + std::to_string(i + 1) + "}", ts.GetDegreeMatrix()[i]);

    MakeTeXFrac(out, "\\lambda_{" + ts.GetName() + "} =  \\bigoplus ^{" + std::to_string(ts.GetMatrix().GetHeight()) +
                "} _{k = 1} \\operatorname{tr} ^{\\frac 1 k} \\left( \\bold{" +
                ts.GetName() + "} \\right)^k", ts.GetSpectral());
    MakeTeXMatrix(out, "\\left(\\lambda_{" + ts.GetName() + "}^{-1}\\bold{" + ts.GetName() + "}\\right)^*", ts.GetKleene());
    MakeTeXMatrix(out, "\\bold{x}_{" + ts.GetName() + "}", ts.GetKleeneWithoutCorrel());
    MakeTeXFrac(out, "\\delta_{" + ts.GetName() + "} = \\bold{1}^T \\left(\\lambda_{" + ts.GetName() +
                "}^{-1}\\bold{" + ts.GetName() + "}\\right)^* \\bold{1}", ts.GetDelta());
    MakeTeXMatrix(out, "\\left(\\delta_{" + ts.GetName() + "}^{-1}\\bold{11}^T \\oplus \\lambda_{" +
                  ts.GetName() + "}^{-1}\\bold{" + ts.GetName() + "}\\right)^*", ts.GetWorstMatrix());

    for (size_t i{}; i < ts.GetWorst().size(); ++i)
    {
        if (ts.GetWorst().size() > 1)
            out << "\t" << i + 1 << ".\n";

        MakeTeXMatrix(out, "\\bold{x}'_{" + ts.GetName() + "}", ts.GetWorst()[i]);
    }

    MakeTeXMatrix(out, "\\bold{P}_{" + ts.GetName() + "}", ts.GetP());

    for (size_t i{}; i < ts.GetPlk().size(); ++i)
    {
        if (ts.GetPlk().size() > 1)
            out << "\t" << i + 1 << ")\n";

        MakeTeXMatrix(out, "\\bold{" + ts.GetName() + "}: \\bold{P}_{lk}", ts.GetPlk()[i]);
        MakeTeXMatrix(out, "\\bold{P}_{" + ts.GetName() + "}\\left(\\bold{I} \\oplus \\bold{P}_{lk}^{-}" +
                        "\\bold{P}_{" + ts.GetName() + "}\\right)", ts.GetBestMatrix()[i]);

        for (size_t j{}; j < ts.GetBest()[i].size(); ++j)
        {
            if (ts.GetBest()[i].size() > 1)
                out << "\t" << j + 1 << ".\n";

            MakeTeXMatrix(out, "\\bold{x}''_{" + ts.GetName() + "}", ts.GetBest()[i][j]);
        }           
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const TropicoSolve& ts)
{
    out << "Matrix " << ts.GetName() << ":\n";
    out << ts.GetMatrix() << "~\n";
    ToDoubleOut(out, ts.GetMatrix()) << "\n";

    for (size_t i{1}; i < ts.GetDegreeMatrix().size(); ++i)
    {
        out << "(" << ts.GetName() << ")^{" << i + 1 << "}:\n" << ts.GetDegreeMatrix()[i];
        ToDoubleOut(out, ts.GetDegreeMatrix()[i]) << "\n";
    }

    out << "Spectral Radius: " << ts.GetSpectral() << " ~ " << ts.GetSpectral().ToDouble() << "\n\n";
    out << "Kleene:\n" << ts.GetKleene() << "~\n";
    ToDoubleOut(out, ts.GetKleene()) << "\n";
    out << "Kleene without correlated cols:\n" << ts.GetKleeneWithoutCorrel() << "~\n";
    ToDoubleOut(out, ts.GetKleeneWithoutCorrel()) << "\n";
    out << "Delta: \n" << ts.GetDelta() << " ~ " << ts.GetDelta().ToDouble() << "\n\n";
    out << "Worst Solve Matrix:\n" << ts.GetWorstMatrix() << "~\n";
    ToDoubleOut(out, ts.GetWorstMatrix()) << "\n";

    for (size_t i{}; i < ts.GetWorst().size(); ++i)
    {
        if (ts.GetWorst().size() > 1)
            out << "\t" << i + 1 << ".\n";

        out << "Worst Solve:\n" << ts.GetWorst()[i] << "~\n";
        ToDoubleOut(out, ts.GetWorst()[i]) << "\n";
    }

    out << "P:\n" << ts.GetP() << "~\n";
    ToDoubleOut(out, ts.GetP()) << "\n";

    for (size_t i{}; i < ts.GetPlk().size(); ++i)
    {
        if (ts.GetPlk().size() > 1)
            out << i + 1 << ")\n";

        out << "P_{lk}:\n" << ts.GetPlk()[i] << "~\n";
        ToDoubleOut(out, ts.GetPlk()[i]) << "\n";
        out << "Best Solve Matrix:\n" << ts.GetBestMatrix()[i] << "~\n";
        ToDoubleOut(out, ts.GetBestMatrix()[i]) << "\n";

        for (size_t j{}; j < ts.GetBest()[i].size(); ++j)
        {
            if (ts.GetBest()[i].size() > 1)
                out << "\t" << j + 1 << ".\n";

            out << "Best Solve:\n" << ts.GetBest()[i][j] << "~\n";
            ToDoubleOut(out, ts.GetBest()[i][j]) << "\n";
        }   
    }

    return out;
}

void TropicoMultiSolve::Solve()
{
    tsC = C;

    tsC.Solve();
    tsC.GiveName("C");

    BWorst.resize(tsC.GetWorst().size());

    for (size_t i{}; i < tsC.GetWorst().size(); ++i)
    {
        BWorst[i].Resize(A[0].GetHeight(), A[0].GetWidth());

        for (size_t j{}; j < tsC.GetWorst()[i].GetHeight(); ++j)
            BWorst[i] += tsC.GetWorst()[i][j][0] * A[j];

        tsBW.push_back(TropicoSolve(BWorst[i]));

        tsBW.back().Solve();

        if (tsC.GetWorst().size() > 1)
            tsBW.back().GiveName("B_1^{(" + std::to_string(i + 1) + ")}");
        else
            tsBW.back().GiveName("B_1");
    }

    BBest.resize(tsC.GetBest().size());

    for (size_t i{}; i < tsC.GetBest().size(); ++i)
    {
        BBest[i].resize(tsC.GetBest()[i].size());

        tsBB.push_back(std::vector<TropicoSolve>());

        for (size_t j{}; j < tsC.GetBest()[i].size(); ++j)
        {
            BBest[i][j].Resize(A[0].GetHeight(), A[0].GetWidth());

            for (size_t t{}; t < tsC.GetBest()[i][j].GetHeight(); ++t)
                BBest[i][j] += tsC.GetBest()[i][j][i][0] * A[t];

            tsBB.back().push_back(TropicoSolve(BBest[i][j]));

            tsBB.back().back().Solve();

            std::string name = "B_2^{";

            if (tsC.GetBest().size() > 1)
                name += "(" + std::to_string(i + 1) + ")";

            if (tsC.GetBest()[i].size() > 1)
                name += "[" + std::to_string(j + 1) + "]";

            name += "}";

            tsBB.back().back().GiveName(name);
        }
    }
}

std::ostream& MakeTeX(std::ostream& out, const TropicoMultiSolve& tms)
{
    MakeTeXMatrix(out, "\\bold{C}", tms.C);

    for (size_t i{}; i < tms.A.size(); ++i)
        MakeTeXMatrix(out, "\\bold{A}_{" + std::to_string(i + 1) + "}", tms.A[i]);

    MakeTeX(out, tms.tsC) << "\n";

    for (size_t i{}; i < tms.tsBW.size(); ++i)
    {
        if (tms.tsBW.size() > 1)
            out << "\t" << i + 1 << "]\n";

        MakeTeXMatrix(out, "\\bold{\\omega}", tms.tsC.GetWorst()[i]) << "\n";
        MakeTeX(out, tms.tsBW[i]) << "\n";
    }

    for (size_t i{}; i < tms.tsBB.size(); ++i)
    {
        if (tms.tsBB.size() > 1)
            out << i + 1 << "]\n";

        for (size_t j{}; j < tms.tsBB[i].size(); ++j)
        {
            if (tms.tsBB[i].size() > 1)
                out << "\t" << j + 1 << ".\n";

            MakeTeXMatrix(out, "", tms.tsC.GetBest()[i][j]) << "\n";
            MakeTeX(out, tms.tsBB[i][j]) << "\n";
        }
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const TropicoMultiSolve& tms)
{
    out << "C:\n" << tms.C << "~\n";
    ToDoubleOut(out, tms.C) << "\n";

    for (size_t i{}; i < tms.A.size(); ++i)
    {
        out << "A_" << i + 1 << "\n" << tms.A[i] << "~\n";
        ToDoubleOut(out, tms.A[i]) << "\n";
    }

    out << tms.tsC << "\n";

    for (size_t i{}; i < tms.tsBW.size(); ++i)
    {
        if (tms.tsBW.size() > 1)
            out << "\t" << i + 1 << ")\n";

        out << tms.tsBW[i] << "\n";
    }

    for (size_t i{}; i < tms.tsBB.size(); ++i)
    {
        if (tms.tsBB.size() > 1)
            out << i + 1 << ")\n";

        for (size_t j{}; j < tms.tsBB[i].size(); ++j)
        {
            if (tms.tsBB[i].size() > 1)
                out << "\t" << j + 1 << ".\n";

            out << tms.tsBB[i][j] << "\n";
        }
    }

    return out;
}

std::istream& operator>>(std::istream& in, TropicoMultiSolve& tms)
{
    TropicoMatrix tm;
    in >> tms.C;

    for (size_t i{}; i < tms.C.GetHeight(); ++i)
    {
        in >> tm;
        tms.A.push_back(tm);
    }

    return in;
}
