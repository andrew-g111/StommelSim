// Author: Andrew Gill
// Very crude testing ground to observe the effects of changing Stommel model's
// parameters. 

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include "Stommel.h"
using namespace std;

void graphEquilibria(int n, double step);

int main() {
    //graphEquilibria(20, 0.05);

    Stommel testStommel = Stommel();
    double testCoords[2] = {1, 2};
    cout << testCoords[0] << "|" << testCoords[1] << endl;
    testStommel.incrementTime(testCoords, 20, 0.01);
    cout << testCoords[0] << "|" << testCoords[1] << endl;
    testStommel.print();
    return 0;
}
    

// Outputs a chart of the number of equilibria, varying over values of d and l.
// n is the number of rows/columns to show, and step is the step size.
void graphEquilibria(int n, double step) {
    Stommel* modelArray[n][n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            modelArray[i][j] = new Stommel((i + 1) * step, 2.0, (j + 1) * step);
        }
    }

    std::cout << "d \\ l";
    for (int i = 0; i < n; i++) {
        std::cout << "|" << std::setw(5) << std::setfill('0') << std::left << (i + 1) * step;
    }
    std::cout << endl;
    for (int i = 0; i < n; i++) {
        std::cout << std::setw(5) << std::setfill('0') << std::left << (i + 1) * step;
        for (int j = 0; j < n; j++) {
            std::cout << "|" << std::setw(5) << std::setfill(' ') << std::right << modelArray[i][j]->numEquilibria();
        }
        std::cout << endl;
    }
}
