#include <iostream>
#include "hw.h"
#include <omp.h>
#include <string.h>

#define N __X_N //10001

int matA[N*N], matB[N*N],
    matCm[N*N], matCm2[N*N];

void calc_4x4(int *result, int *A, int *B, int len)
{
    //B转置过
    register int c00_reg=0, c01_reg=0, c02_reg=0, c03_reg=0;
    register int c10_reg=0, c11_reg=0, c12_reg=0, c13_reg=0;
    register int c20_reg=0, c21_reg=0, c22_reg=0, c23_reg=0;
    register int c30_reg=0, c31_reg=0, c32_reg=0, c33_reg=0;

    int *b0_ptr=&B[0];
    int *b1_ptr=&B[len];
    int *b2_ptr=&B[2*len];
    int *b3_ptr=&B[3*len];

    int *a0_ptr=&A[0];
    int *a1_ptr=&A[len];
    int *a2_ptr=&A[2*len];
    int *a3_ptr=&A[3*len];

    register int a0_reg, a1_reg, a2_reg, a3_reg;
    register int b0_reg, b1_reg, b2_reg, b3_reg;

    for(int i=0; i<len; i++)
    {
        a0_reg=*a0_ptr++;
        a1_reg=*a1_ptr++;
        a2_reg=*a2_ptr++;
        a3_reg=*a3_ptr++;

        b0_reg=*b0_ptr++;
        b1_reg=*b1_ptr++;
        b2_reg=*b2_ptr++;
        b3_reg=*b3_ptr++;

        c00_reg+=a0_reg*b0_reg;
        c01_reg+=a0_reg*b1_reg;
        c02_reg+=a0_reg*b2_reg;
        c03_reg+=a0_reg*b3_reg;

        c10_reg+=a1_reg*b0_reg;
        c11_reg+=a1_reg*b1_reg;
        c12_reg+=a1_reg*b2_reg;
        c13_reg+=a1_reg*b3_reg;

        c20_reg+=a2_reg*b0_reg;
        c21_reg+=a2_reg*b1_reg;
        c22_reg+=a2_reg*b2_reg;
        c23_reg+=a2_reg*b3_reg;

        c30_reg+=a3_reg*b0_reg;
        c31_reg+=a3_reg*b1_reg;
        c32_reg+=a3_reg*b2_reg;
        c33_reg+=a3_reg*b3_reg;
    }

    result[0]=c00_reg;
    result[1]=c01_reg;
    result[2]=c02_reg;
    result[3]=c03_reg;

    result[len+0]=c10_reg;
    result[len+1]=c11_reg;
    result[len+2]=c12_reg;
    result[len+3]=c13_reg;

    result[2*len+0]=c20_reg;
    result[2*len+1]=c21_reg;
    result[2*len+2]=c22_reg;
    result[2*len+3]=c23_reg;

    result[3*len+0]=c30_reg;
    result[3*len+1]=c31_reg;
    result[3*len+2]=c32_reg;
    result[3*len+3]=c33_reg;

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
        for (int i=0; i<N; i+=4)
            for (int j=0; j<N; j+=4)
                calc_4x4(&matC2[i*N+j], &matC[i*N], &matA[j*N], N);

        int *t = matC; matC = matC2; matC2 = t;
    }

    cout<<"calc over"<<endl;
    output(matC, n);
}
