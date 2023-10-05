#include <iostream>
#include "hw.h"
#include <omp.h>
#include <string.h>

#define N __X_N //10001

int matA[N*N], matB[N*N],
    matCm[N*N], matCm2[N*N];

void calc_1x4(int *result, int *A, int *B, int len)
{
    //B转置过
    result[0]=0;
    result[1]=0;
    result[2]=0;
    result[3]=0;
    // memset(result, 0, sizeof(int[4]));

    for(int i=0; i<len; i++)
    {
        result[0]+=A[i]*B[i];
        result[1]+=A[i]*B[len+i];
        result[2]+=A[i]*B[2*len+i];
        result[3]+=A[i]*B[3*len+i];
    }
}

int main()
{
    using namespace std;

    input(matA, matB);  //随机生成A和B
    int *matC = matCm, *matC2 = matCm2, n = 1;
    memcpy(matC, matA, sizeof(int[N*N]));

    //转置A
    #pragma omp parallel for schedule(guided, 100)
    for (int i=0;i<N;++i)
        for (int j=0;j<i;++j) {
            int t = matA[i*N+j];
            matA[i*N+j] = matA[j*N+i];
            matA[j*N+i] = t;
        }

    //转置B
    #pragma omp parallel for schedule(guided, 100)
    for (int i=0;i<N;++i)
        for (int j=0;j<i;++j) {
            int t = matB[i*N+j];
            matB[i*N+j] = matB[j*N+i];
            matB[j*N+i] = t;
        }

    for (int k=0;k<n;++k)
    {
        cout<<"n: "<<n<<"; k: "<<k<<endl;

        #pragma omp parallel for
        for (int i=0;i<N*N;++i)
            matA[i] += matB[i];
        
        // memset(matC2, 0, sizeof(int[N*N]));

        #pragma omp parallel for
        for (int i=0; i<N; ++i)
            for (int j=0; j<N; j+=4)
            {
                //计算第i行第j列的结果
                calc_1x4(&matC2[i*N+j], &matC[i*N], &matA[j*N], N);
            }

        int *t = matC; matC = matC2; matC2 = t;
    }

    cout<<"calc over"<<endl;
    output(matC, n);
}
