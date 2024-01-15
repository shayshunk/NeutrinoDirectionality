#include <iostream>

std::ostream& boldOn(std::ostream& os)
{
    return os << "\033[1m";
}

std::ostream& underlineOn(std::ostream& os)
{
    return os << "\033[4m";
}

std::ostream& redOn(std::ostream& os)
{
    return os << "\033[31m";
}

std::ostream& greenOn(std::ostream& os)
{
    return os << "\033[32m";
}

std::ostream& yellowOn(std::ostream& os)
{
    return os << "\033[33m";
}

std::ostream& blueOn(std::ostream& os)
{
    return os << "\033[34m";
}

std::ostream& cyanOn(std::ostream& os)
{
    return os << "\033[36m";
}

std::ostream& whiteOn(std::ostream& os)
{
    return os << "\033[37m";
}

std::ostream& resetFormats(std::ostream& os)
{
    return os << "\033[0m";
}