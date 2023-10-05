#include <iostream>
#include "hw.h"
#include <omp.h>

#define N __X_N //10001

int matA[N*N], matB[N*N],
    matCm[N*N], matCm2[N*N];

int main()
{
    using namespace std;

    input(matA, matB);  //随机生成A和B
    int *matC = matCm, *matC2 = matCm2, n = 1;
    memcpy(matC, matA, sizeof(int[N*N]));

// //转置A
// //#pragma omp parallel for schedule(guided, 100)
//     for (int i=0;i<N;++i)
//         for (int j=0;j<i;++j) {
//             int t = matA[i*N+j];
//             matA[i*N+j] = matA[j*N+i];
//             matA[j*N+i] = t;
//         }

// //转置B
// //#pragma omp parallel for schedule(guided, 100)
//     for (int i=0;i<N;++i)
//         for (int j=0;j<i;++j) {
//             int t = matB[i*N+j];
//             matB[i*N+j] = matB[j*N+i];
//             matB[j*N+i] = t;
//         }
    
    for (int k=0;k<n;++k)
    {
        cout<<"n: "<<n<<"; k: "<<k<<endl;
        #pragma omp parallel for
        for (int i=0;i<N*N;++i)
            matA[i] += matB[i];

        #pragma omp parallel for
        for (int i=0;i<N;++i)
        {
            // cout<<i<<endl;

            for (int j=0;j<N;++j)
            {
                int sum = 0;
                for (int k=0;k<N;++k)
                    sum += matC[i*N+k] * matA[k*N+j];
                matC2[i*N+j] = sum;
            }
        }

        int *t = matC; matC = matC2; matC2 = t;
    }

    cout<<"calc over"<<endl;
    output(matC, n);
}
