#include <iostream>
#include <iomanip>
#include <cmath>
#include "InverseCumulativeNormal.h"

int main() {
    quant::InverseCumulativeNormal icn;

    double xs[] = {1e-12, 1e-6, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999999, 1.0};
    for (double x : xs) {
        double z = icn(x);
        std::cout << std::setprecision(17)
                  << "x =" << std::setw(18) << std::defaultfloat << x
                  << "  (hex: " << std::hexfloat << x << std::defaultfloat << ")"
                  << "  z =" << std::setw(14) << std::setprecision(10) << z
                  << "\n";
    }
    return 0;
}
