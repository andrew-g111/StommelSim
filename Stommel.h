// Author: Andrew Gill
// Defines the Stommel class, which is initialized with lambda, R, and delta.
// A Stommel object can generate equilibrium points and test their stability.

#include<cmath>
#include<vector>
#include<complex>
#include<iostream>
#include<iomanip>
using namespace std;

class Stommel {

    public:
        // Default values: delta is 1/6, r is 2, and lambda is 1/5.
        Stommel();

        // Initializes constants for an instance of the Stommel model. As
        // such, the class may not work correctly if inputD <= 0.
        Stommel(double inputD, double inputR, double inputL);

        // Getter for equilibria
        vector<double> getEquilibria() {return equilibria;}

        // Returns the number of equilibria.
        int numEquilibria() {return equilibria.size();}

        // Prints lambda, r, delta, and the contents of equilibria.
        void print();

        // Updates the array coords = {x, y} with how the system evolves after
        // numIncrements * dt time. dt controls how large each step is, while
        // numIncrements is the number of steps.
        void incrementTime(double coords[2], int numIncrements, double dt);

        // Returns the flow at x, y.
        double flowAt(double x, double y) {return (-y + r * x) / lambda;}

        // Returns the flow at coords = {x, y}
        double flowAt(double coords[2]) {return flowAt(coords[0], coords[1]);}

        // Returns an array {x, y} of the system's state at the equilibrium
        // associated with the given f. Assumes f is a valid equilibrium.
        double* coordsAtEquilibrium(double f);
        
    protected:
        double delta, 
               r, 
               lambda;

        // equilibrium values of f
        vector<double> equilibria;

        // Returns a vector containing up to six candidates for equilibria.
        // Notably, this vector contains every equilibrium point.
        vector<double> getCandidates();

        // Checks each element of candidates, returning a vector with only
        // the elements of candidates that are equilibria.
        vector<double> trimCandidates(vector<double> candidates);

        // Simple bubble sort in place for the list of equilibria. Not to be
        // used for any large lists, but it's impossible to have more than
        // six equilibria, so it doesn't matter that much here.
        void bubbleSort(vector<double> &list);

    private:
        const double TINY_NUM = .00001; // Margin of error for floating point
                                        // math, shouldn't affect accuracy of
                                        // equilibrium values.

        

         // Returns the difference between the two sides of the condition
        // for equilibrium. The closer f is to equilibrium, the closer to
        // zero this function will return.
        double equilibriumEpsilon(double f);

        // Returns a vector containing the real roots of a cubic polynomial
        // ax^3+bx^2+cx+d where a is nonzero.
        vector<double> realRoots(double a, double b, double c, double d);

        
        
};


Stommel::Stommel() {
    delta = 1.0 / 6.0; 
    r = 2.0; 
    lambda = 1.0 / 5.0;
    equilibria = trimCandidates(getCandidates());
    bubbleSort(equilibria);
}

   
Stommel::Stommel(double inputD, double inputR, double inputL) {
    delta = inputD;
    r = inputR;
    lambda = inputL;
    equilibria = trimCandidates(getCandidates());
    bubbleSort(equilibria);
}


double Stommel::equilibriumEpsilon(double f) {
    double frac1 = -1.0 / (1.0 + fabs(f)),
        frac2 = r / (1.0 + (fabs(f) / delta));
    return fabs((lambda * f) - (frac1 + frac2));
}


vector<double> Stommel::realRoots(double a, double b, double c, double d) {
    vector<double> realRoots;
    double delta0 = (b * b) - (3.0 * a * c);
    double delta1 = (2 * pow(b, 3)) - (9 * a * b * c) + (27 * a * a * d);
    complex<double> deltaExpression = pow(delta1, 2) - 4 * pow(delta0, 3);
    complex<double> cRadicand = (delta1 + sqrt(deltaExpression)) / 2.0;

    complex<double> bigC = pow(cRadicand, 1.0 / 3.0);

    complex<double> root = -1.0 * (b + bigC + (delta0 / bigC)) / (3 * a);
    complex<double> rootOfUnity = {-1.0 / 2.0, sqrt(3.0) / 2.0};

    for (int i = 0; i < 3; i++) {   // Search all 3 roots
        if (fabs(root.imag()) < TINY_NUM) {
            realRoots.push_back(root.real());
        }
        bigC *= rootOfUnity;
        root = -1.0 * (b + bigC + (delta0 / bigC)) / (3 * a);
    }

    if (!realRoots.size()) {
        throw logic_error("realRoots should return at least one root");
    }
    
    return realRoots;

}


vector<double> Stommel::getCandidates() {
    vector<double> candidates;
    vector<double> roots1 = realRoots(lambda, 
                                     lambda * (1 + delta),
                                     delta * (lambda - r) + 1,
                                     delta * (1 - r));
    vector<double> roots2 = realRoots(lambda, 
                                     -lambda * (1 + delta),
                                     delta * (lambda + r) - 1,
                                     delta * (1 - r));

    for (int i = 0; i < roots2.size(); i++) {
        roots1.push_back(roots2.at(i));
    }

    return roots1;
}


vector<double> Stommel::trimCandidates(vector<double> candidates) {
    vector<double> finalEquilibria;
    for (int i = 0; i < candidates.size(); i++) {
        if (equilibriumEpsilon(candidates[i]) < TINY_NUM) {
            finalEquilibria.push_back(candidates[i]);
        }
    }

    return finalEquilibria;
}


void Stommel::print() {
    cout << setprecision(10)
         << "lambda: " << lambda << endl
         << "R: " << r << endl
         << "delta: " << delta << endl
         << "equilibrium values of f:" << endl;
    
    for (int i = 0; i < equilibria.size(); i++) {
        cout << equilibria.at(i) << endl;
    }
}


void Stommel::incrementTime(double coords[2], int numIncrements, double dt) {
    
    double x = coords[0],
           y = coords[1],
           dxdt,
           dydt;

    for (int i = 0; i < numIncrements; i++) {
       dxdt = delta * (1.0 - x) - (x / lambda) * fabs(-y + r * x);
       dydt = 1.0 - y - (y / lambda) * fabs(-y + r * x);
       x += dxdt * dt;
       y += dydt * dt;
    }

    coords[0] = x;
    coords[1] = y;
}


double* Stommel::coordsAtEquilibrium(double f) {
    static double coords[2];
    coords[0] = 1.0 / (1.0 + fabs(f) / delta);
    coords[1] = 1.0 / (1.0 + fabs(f));
    return coords;
}


void Stommel::bubbleSort(vector<double> &list) {
    int size = list.size();
    double temp;
    if (size < 2) { // Prevents indexing errors
        return;
    }

    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - 1 - i; j++) {
            if (list[j] > list[j + 1]) {
                temp = list[j];
                list[j] = list[j+1];
                list[j+1] = temp;
            }
        }
    }
}