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
    register int a_reg;
    register int c0_reg, c1_reg, c2_reg, c3_reg, c4_reg, c5_reg, c6_reg, c7_reg;
    c0_reg=0;
    c1_reg=0;
    c2_reg=0;
    c3_reg=0;
    c4_reg=0;
    c5_reg=0;
    c6_reg=0;
    c7_reg=0;

    int *b0_ptr=&B[0];
    int *b1_ptr=&B[len];
    int *b2_ptr=&B[2*len];
    int *b3_ptr=&B[3*len];
    int *b4_ptr=&B[4*len];
    int *b5_ptr=&B[5*len];
    int *b6_ptr=&B[6*len];
    int *b7_ptr=&B[7*len];

    for(int i=0; i<len; i++)
    {
        a_reg=A[i];

        c0_reg+=a_reg* *b0_ptr++;
        c1_reg+=a_reg* *b1_ptr++;
        c2_reg+=a_reg* *b2_ptr++;
        c3_reg+=a_reg* *b3_ptr++;
        c4_reg+=a_reg* *b4_ptr++;
        c5_reg+=a_reg* *b5_ptr++;
        c6_reg+=a_reg* *b6_ptr++;
        c7_reg+=a_reg* *b7_ptr++;
    }

    result[0]=c0_reg;
    result[1]=c1_reg;
    result[2]=c2_reg;
    result[3]=c3_reg;
    result[4]=c4_reg;
    result[5]=c5_reg;
    result[6]=c6_reg;
    result[7]=c7_reg;

    return;
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
            for (int j=0; j<N; j+=8)
            {
                //计算第i行第j列的结果
                calc_1x4(&matC2[i*N+j], &matC[i*N], &matA[j*N], N);
            }

        int *t = matC; matC = matC2; matC2 = t;
    }

    cout<<"calc over"<<endl;
    output(matC, n);
}
