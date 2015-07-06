

#include <iostream>

using namespace std;

#define NMAX 21

struct polinom
{
    double coefs[NMAX];
    int pow;
};





int main() {

    int N;
    cin >> N;

    polinom p;
    p.pow = N;

    for (int i = 0; i <= N; ++i) {
        cin >> p.coefs[i];
    }

    int sum = (p.coefs[p.pow-1] / p.coefs[p.pow]) * -1;
    int product = (p.coefs[0] / p.coefs[p.pow]) * pow(-1., p.pow);


    cout << sum << " " << product << endl;

    return 0;
}
