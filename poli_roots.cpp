

#include <iostream>
#include <cmath>

using namespace std;


const int NMAX = 21;

typedef int INT;
typedef INT (COEFS) [NMAX];
typedef struct {
    COEFS coefs;
    int pow;
} POLINOM;


ostream& operator<< (ostream& os, const POLINOM& p) {
    cout << "polinom: degree " << p.pow << ", coefs {";

    for (int i = 0; i <= p.pow; ++i)
        cout << p.coefs[i] << ",";
    cout << "}";

    return os;
}

double calc_poli(POLINOM& p, double x) {
    double res = p.coefs[0];
    double power = x;

    if (1 <= p.pow) {
        res += p.coefs[1] * x;
        power *= x;
    

        for (int i = 2; i <= p.pow; ++i) {
            res += p.coefs[i] * power;
            power *= x;
        }
    }
    
    return res;
}


void solve(POLINOM& p, int& sum, int& product) {
    const double e = 0.1;

    POLINOM dp;
    if (2 < p.pow) {
        dp.pow = p.pow - 1;
        for (int i = 1; i <= p.pow; ++i)
            dp.coefs[i-1] = p.coefs[i] * i;
    }

    int root_number = p.pow;

    for (int r = 0; r < root_number; ++r) {
        if (1 == p.pow) {
            int root = -p.coefs[0] / p.coefs[1];
            sum += root;
            product *= root;
            break;
        }
        else if (2 == p.pow) {
            double d = sqrt(p.coefs[1] * p.coefs[1] - 4. * p.coefs[0] * p.coefs[2]);
            int root1 = (-p.coefs[1] + d) / (p.coefs[2] << 1);
            int root2 = (-p.coefs[1] - d) / (p.coefs[2] << 1);
            sum += root1 + root2;
            product *= root1 * root2;
            break;
        }
        else {
            double x0 = 400.;
            double x1;
            double res = calc_poli(p, x0);

            while ((res < 0 ? -res : res) > e) {
                x1 = x0 - res / calc_poli(dp, x0);
                res = calc_poli(p, x1);
                x0 = x1;
            }
          
            int root = (int)round(x1); 
            sum += root;
            product *= root;

            p.coefs[0] /= -root; 
            for (int i = 1; i < p.pow; ++i) {
                p.coefs[i] = (p.coefs[i] - p.coefs[i-1]) / (-root);
            }
            --p.pow;

            dp.pow = p.pow - 1;
            for (int i = 1; i <= p.pow; ++i)
                dp.coefs[i-1] = p.coefs[i] * i;
        }
    }
}


int main(int argc, const char* argv[]) {

    int N = 0;
    cin >> N;

    POLINOM p;
    p.pow = N;
    for (int i = 0; i <= p.pow; ++i)
        cin >> p.coefs[i];

    int sum = 0;
    int product = 1;
    solve(p, sum, product);

    cout << sum << " " << product << endl;

    return 0;
}

