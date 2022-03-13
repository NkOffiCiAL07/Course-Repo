#include <iostream> // Including necessary libraries.
#include <tuple>
using namespace std;
int steps = 0;

/*Function Description-
    input: int a, int b
    output: tuple<int, int, int> i.e. returns (x, y, gcd) respectively
    info: This function calculates gcd of given two numbers and coefficients x,y such that a*x + b*y = gcd(a,b)
*/
tuple<int, int, int> gcd(int a, int b)
{
    steps++;
    if (a == 0) // base step.
        return make_tuple(0, 1, b);

    tuple<int, int, int> value = gcd(b % a, a);                       // coefficient calculation using recursive step.
    return make_tuple((get<1>(value) - (b / a) * get<0>(value)), get<0>(value), get<2>(value)); // returning as a tuple.
}

// driver function code.
int main()
{
    int a, b;
    cout << "Input the two numbers : ";
    cin >> a >> b;

    tuple<int, int, int> coefficient = gcd(a, b); // calling extented_coefficient function using tuple
    int coef1 = get<0>(coefficient);
    int coef2 = get<1>(coefficient);
    int gcd = get<2>(coefficient);

    cout << "The GCD of " << a << " and " << b << " is " << gcd << endl; // Outputs.
    cout << "x = " << coef1 << ","
         << " y = " << coef2 << endl;
    cout << gcd << " = " << a << "*" << coef1 << " + " << b << "*" << coef2 << endl; 
    cout<<"steps = "<<steps<<"\n";

    return 0;
}
