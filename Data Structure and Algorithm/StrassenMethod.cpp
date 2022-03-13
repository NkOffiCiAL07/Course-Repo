// In This code we randomly generate two n*n matrices (size n = 2^k) compare it with n^log7(or 7^k)
// and output the running time for n, n^log 7, and the ratio of both these values for n

#include <iostream> // including libraries
#include <vector>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono;

// subtracting matrices
void matSubtr(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C, int sz)
{
    for (int i = 0; i < sz; i++)
    {
        for (int j = 0; j < sz; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}
// adding matrices
void matAdd(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C, int sz)
{
    for (int i = 0; i < sz; i++)
    {
        for (int j = 0; j < sz; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// recursive multiplication function using Strassen's method
void strassMltp(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C, int sz)
{
    if (sz == 1) // base condition
    {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }
    int new_size = sz / 2; // new size of matrix reduced to half

    vector<int> zeroVec(new_size, 0); // creating zero matrix

    // creating new matrices and initialize these with zer0 matrix
    vector<vector<int>> A11(new_size, zeroVec);
    vector<vector<int>> A12(new_size, zeroVec);
    vector<vector<int>> A21(new_size, zeroVec);
    vector<vector<int>> A22(new_size, zeroVec);

    vector<vector<int>> B11(new_size, zeroVec);
    vector<vector<int>> B12(new_size, zeroVec);
    vector<vector<int>> B21(new_size, zeroVec);
    vector<vector<int>> B22(new_size, zeroVec);

    vector<vector<int>> C11(new_size, zeroVec);
    vector<vector<int>> C12(new_size, zeroVec);
    vector<vector<int>> C21(new_size, zeroVec);
    vector<vector<int>> C22(new_size, zeroVec);

    vector<vector<int>> temp1(new_size, zeroVec);
    vector<vector<int>> temp2(new_size, zeroVec);
    vector<vector<int>> temp3(new_size, zeroVec);
    vector<vector<int>> temp4(new_size, zeroVec);
    vector<vector<int>> temp5(new_size, zeroVec);
    vector<vector<int>> temp6(new_size, zeroVec);
    vector<vector<int>> temp7(new_size, zeroVec);
    vector<vector<int>> temp8(new_size, zeroVec);
    vector<vector<int>> temp9(new_size, zeroVec);
    vector<vector<int>> temp10(new_size, zeroVec);

    vector<vector<int>> P1(new_size, zeroVec);
    vector<vector<int>> P2(new_size, zeroVec);
    vector<vector<int>> P3(new_size, zeroVec);
    vector<vector<int>> P4(new_size, zeroVec);
    vector<vector<int>> P5(new_size, zeroVec);
    vector<vector<int>> P6(new_size, zeroVec);
    vector<vector<int>> P7(new_size, zeroVec);

    // temporary matrix for storing values
    vector<vector<int>> mA(new_size, zeroVec);
    vector<vector<int>> mB(new_size, zeroVec);

    // filling values of A11,A12,....
    // or segmenting the given A and B matrices
    for (int i = 0; i < new_size; i++)
    {
        for (int j = 0; j < new_size; j++)
        {
            A11[i][j] = A[i][j];
            B11[i][j] = B[i][j];

            A12[i][j] = A[i][j + new_size];
            B12[i][j] = B[i][j + new_size];

            A21[i][j] = A[i + new_size][j];
            B21[i][j] = B[i + new_size][j];

            A22[i][j] = A[i + new_size][j + new_size];
            B22[i][j] = B[i + new_size][j + new_size];
        }
    }

    matSubtr(B12, B22, temp1, new_size); // temp1 = B12 - B22
    matAdd(A11, A12, temp2, new_size);   // temp2 = A11 + A12
    matAdd(A21, A22, temp3, new_size);   // temp3 = A21 + A22
    matSubtr(B21, B11, temp4, new_size); // temp4 = B21 - B11
    matAdd(A11, A22, temp5, new_size);   // temp5 = A11 + A22
    matAdd(B11, B22, temp6, new_size);   // temp6 = B11 + B22
    matSubtr(A12, A22, temp7, new_size); // temp7 = A12 - A22
    matAdd(B21, B22, temp8, new_size);   // temp8 = B21 + B22
    matSubtr(A11, A21, temp9, new_size); // temp9 = A11 - A21
    matAdd(B11, B12, temp10, new_size);  // temp10 = B11 + B12

    strassMltp(A11, temp1, P1, new_size);    // P1 = A11 * temp1
    strassMltp(temp2, B22, P2, new_size);    // P2 = temp2 * B22
    strassMltp(temp3, B11, P3, new_size);    // P3 = temp3 * B11
    strassMltp(A22, temp4, P4, new_size);    // P4 = A22 * temp4
    strassMltp(temp5, temp6, P5, new_size);  // P5 = temp5 * temp6
    strassMltp(temp7, temp8, P6, new_size);  // P6 = temp7 * temp8
    strassMltp(temp9, temp10, P7, new_size); // P7 = temp9 * temp10

    // C11 = P5 + P4 - P2 + P6
    matAdd(P5, P4, mA, new_size);    // mA = P5 + P4
    matAdd(mA, P6, mB, new_size);    // mB = (P5 + P4) + P6 (or mA + P6)
    matSubtr(mB, P2, C11, new_size); // C11 = (P5 + P4 + P6) - P2
    matAdd(P1, P2, C12, new_size);   // C12 = P1 + P2
    matAdd(P3, P4, C21, new_size);   // C21 = P3 + P4

    // C22 = P5 + P1 - P3 + P7
    matAdd(P5, P1, mA, new_size);    // P5 + P1
    matSubtr(mA, P3, mB, new_size);  // (P5 + P1) - P3
    matSubtr(mB, P7, C22, new_size); // C22 = (P5 + P1 - P3) - P7

    // filling values in C
    for (int i = 0; i < new_size; i++)
    {
        for (int j = 0; j < new_size; j++)
        {
            C[i][j] = C11[i][j];
            C[i + new_size][j] = C21[i][j];
            C[i][j + new_size] = C12[i][j];
            C[i + new_size][j + new_size] = C22[i][j];
        }
    }
}

// printing the given matrix
void printMatrix(vector<vector<int>> &mat, int sz)
{
    for (int i = 0; i < sz; i++)
    {
        for (int j = 0; j < sz; j++)
        {
            cout << mat[i][j] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}
// This function generates matrices A, B with random (0, 1) values and C with zeros
void generate(vector<vector<int>> &A, vector<vector<int>> &B, vector<vector<int>> &C, int n)
{
    // generating matrix A
    for (int i = 0; i < n; ++i)
    {
        vector<int> temp;
        for (int j = 0; j < n; ++j)
        {
            temp.push_back(rand() % 2);
        }
        A.push_back(temp);
    }
    // generating matrix B
    for (int i = 0; i < n; ++i)
    {
        vector<int> temp;
        for (int j = 0; j < n; ++j)
        {
            temp.push_back(rand() % 2);
        }
        B.push_back(temp);
    }
    // generating matrix C
    for (int i = 0; i < n; ++i)
    {
        vector<int> temp;
        for (int j = 0; j < n; ++j)
        {
            temp.push_back(0);
        }
        C.push_back(temp);
    }
}
int main()
{
    int n, k = 1;
    cout<<"Value of n\tTime Taken\tn^log7\t\tRatio (Time Taken/n^log7)\n";
    cout<<"\t\t(in microsec)\n";
    for (k = 1; k <= 8; k++)
    {
        n = pow(2,k);   // n = 2^k
        vector<vector<int>> A, B, C;
        generate(A, B, C, n);   // generating matrices A, B with random (0, 1) values and C with zeros

        auto start = high_resolution_clock::now();  //time calculation
        strassMltp(A, B, C, n);                     // calling function
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout<<n<<"\t\t";
        cout << duration.count() << "\t\t";
        cout << pow(7,k)<<"\t\t";
        cout << (double)duration.count()/pow(7,k)<<"\n";
    }

    return 0;
}
