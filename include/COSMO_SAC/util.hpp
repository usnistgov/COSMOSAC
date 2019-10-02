#ifndef COSMO_UTIL_HPP
#define COSMO_UTIL_HPP

#include <fstream>
#include <sstream>
#include <locale>
#include <vector>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace COSMOSAC{

struct delim : std::numpunct<char> {
    char m_c;
    delim(char c) : m_c(c) {};
    char do_decimal_point() const { return m_c; }
};

/** From CoolProp...
*/
static double mystrtod(const std::string &s, const char c = '.') {
    std::stringstream ss(s); double f;
    ss.imbue(std::locale(ss.getloc(), new delim(c)));
    //ss.setf(std::ios::uppercase);
    ss >> f;
    auto chars_remaining = ss.rdbuf()->in_avail();
    /*if (chars_remaining != 0) {
    throw std::invalid_argument("string [" + s + "] was not converted fully to a float with the punct [" + std::string(1, c) + "]");
    }*/
    return f;
}

/** Based on mystrtod, taken from CoolProp
*/
static int mystrtoi(const std::string &s) {
    std::stringstream ss(s); int i;
    ss >> i;
    if (ss.rdbuf()->in_avail() != 0) {
        throw std::invalid_argument("string [" + s + "] was not converted fully to an int");
    }
    return i;
}

/**  Read in an entire file in one shot
*/
static inline std::string get_file_contents(const std::string &filename) {
    using std::ios;
    std::ifstream ifs(filename.c_str(), ios::in | ios::binary);
    if (!ifs) {
        throw std::invalid_argument("filename cannot be opened: " + filename);
    }
    // See https://stackoverflow.com/a/116220/1360263
    return static_cast<std::stringstream const&>(std::stringstream() << ifs.rdbuf()).str();
}

/** Take in a string, use the delimiter to split it up into parts
*  Inspired by http://www.cplusplus.com/faq/sequences/strings/split/#getline
*/
static std::vector<std::string> str_split(const std::string &s,
    const std::string & delimiters = "\n") {
    int current;
    int next = -1;
    std::vector<std::string> o;
    do {
        current = next + 1;
        next = static_cast<int>(s.find_first_of(delimiters, current));
        o.push_back(s.substr(current, next - current));
    } while (next != std::string::npos);
    return o;
}

/// The following code for the trim functions was taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
    // trim from start
static std::string& strlstrip(std::string& s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}
// trim from end
static std::string& strrstrip(std::string& s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}
// trim from both ends
static std::string& strstrip(std::string& s) {
    return strlstrip(strrstrip(s));
}

} /* namespace COSMOSAC */
#endif
