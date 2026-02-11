#ifndef UTILITY_H
#define UTILITY_H
#include <sstream>
#include <vector>
#include <string>
#include <algorithm> 
#include "my_random.h"

inline double clamp01(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int to_int(double x) { return int(pow(clamp01(x), 1 / 2.2) * 255 + 0.5); }

static std::vector<std::string> split(const std::string& str) {
    std::vector<std::string> result;
    std::string cleaned = str;

    cleaned.erase(std::remove(cleaned.begin(), cleaned.end(), '\r'), cleaned.end());
    cleaned.erase(std::remove(cleaned.begin(), cleaned.end(), '\n'), cleaned.end());
    
    std::istringstream iss(cleaned);
    std::string token;
    
    while (iss >> token)
        result.push_back(token);
    
    return result;
}

inline double clamp(double val, double low, double high) {
	if (val < low)
		return low;
	else if (val > high)
		return high;
	else
		return val;
}

inline double diff_of_products(double a, double b, double c, double d) {
	double cd = c * d;
	double diff_of_prod = fma(a, b, -cd);
	double err = fma(-c, d, cd);
	return diff_of_prod + err;
}

inline double lerp(double x, double a, double b) {
	return (1 - x) * a + x * b;
}

// Gaussian Random Number Generation (Box-Muller transform)
void box_muller(double u1, double u2, double& z0, double& z1) {
	double r = sqrt(-2.0 * log(u1));
	double theta = 2.0 * M_PI * u2;
	z0 = r * cos(theta);
	z1 = r * sin(theta);
}

#endif
