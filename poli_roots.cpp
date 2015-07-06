

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <algorithm>

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

double calc_poli(const POLINOM& p, double x) {
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


double find_initial(const POLINOM& p) {
    srand(time(NULL));

    double left = 3.1415926;
    double result = calc_poli(p, left);

    double beg  = -1000000000.;
    double range = 10000000.;
    
    double right = ((double)rand() / RAND_MAX) * range - beg;
    beg += range;
    double tmp = calc_poli(p, right);

    double best = left;
    double best_val = result;
    if (fabs(result) > fabs(tmp)) {
        best = right;
        best_val = tmp;
    }
   
    const int max_count = 100; 
    int count = 0;
    while ((result > 0) == (tmp > 0) && count++ < max_count) {
        right = ((double)rand() / RAND_MAX) * range - beg;
        beg = beg + range;
        beg = beg > (1000000000 - range) ? -1000000000 : beg;
        tmp = calc_poli(p, right);

        if (fabs(tmp) < fabs(best_val)) {
            best = right;
            best_val = tmp;
        }
    }

    if (count >= max_count)
        return best;

    bool left_positive = result > 0;
    bool right_positive = tmp > 0;


    //
    if (left > right) {
       swap(left, right);
       swap(left_positive, right_positive);
    }

    const double e = 1;

    double mid;
    while ((tmp > 0 ? tmp : -tmp) > e) {
        mid = (right + left) / 2.;
        tmp = calc_poli(p, mid);


        if (left_positive == (tmp > 0))
            left = mid;
        else 
            right = mid;
        
    }

    return mid;
}


void solve(POLINOM& p, int& sum, int& product) {
    const double e = 0.05;

    POLINOM dp;
    if (2 < p.pow) {
        dp.pow = p.pow - 1;
        for (int i = 1; i <= p.pow; ++i)
            dp.coefs[i-1] = p.coefs[i] * i;
    }
//int cnt = 0;
//vector<int> roots;
    int root_number = p.pow;

    for (int r = 0; r < root_number; ++r) {
        if (1 == p.pow) {
            int root = -p.coefs[0] / p.coefs[1];
//cout << "Root " << cnt++ << ": " << root << endl;
//roots.push_back(root);
            sum += root;
            product *= root;
            break;
        }
        else if (2 == p.pow) {
            double d = sqrt(p.coefs[1] * p.coefs[1] - 4. * p.coefs[0] * p.coefs[2]);
            int root1 = (-p.coefs[1] + d) / (p.coefs[2] << 1);
            int root2 = (-p.coefs[1] - d) / (p.coefs[2] << 1);
//cout << "Root^2 " << (cnt++) << ": " << root1 << endl
//     << "Root^2 " << (cnt++) << ": " << root2 << endl;
//roots.push_back(root1);
//roots.push_back(root2);
            sum += root1 + root2;
            product *= root1 * root2;
            break;
        }
        else {
            double x0 = find_initial(p);
            double x1 = x0;
            double res = calc_poli(p, x0);

            while ((res < 0 ? -res : res) > e) {
                x1 = x0 - res / calc_poli(dp, x0);
                res = calc_poli(p, x1);
                x0 = x1;
            }
          
            int root = (int)round(x1); 
//cout << "Root " << cnt++ << ": " << root << endl;
//roots.push_back(root);
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

//sort(roots.begin(), roots.end());
//for (int i = 0; i < roots.size(); ++i)
//    cout << roots[i] << ", ";
//cout << endl;
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

