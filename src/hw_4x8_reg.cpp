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
    register int a0_reg, a1_reg, a2_reg, a3_reg;

    register int c00_reg=0, c01_reg=0, c02_reg=0, c03_reg=0, c04_reg=0, c05_reg=0, c06_reg=0, c07_reg=0;
    register int c10_reg=0, c11_reg=0, c12_reg=0, c13_reg=0, c14_reg=0, c15_reg=0, c16_reg=0, c17_reg=0;
    register int c20_reg=0, c21_reg=0, c22_reg=0, c23_reg=0, c24_reg=0, c25_reg=0, c26_reg=0, c27_reg=0;
    register int c30_reg=0, c31_reg=0, c32_reg=0, c33_reg=0, c34_reg=0, c35_reg=0, c36_reg=0, c37_reg=0;

    int *b0_ptr=&B[0];
    int *b1_ptr=&B[len];
    int *b2_ptr=&B[2*len];
    int *b3_ptr=&B[3*len];
    int *b4_ptr=&B[4*len];
    int *b5_ptr=&B[5*len];
    int *b6_ptr=&B[6*len];
    int *b7_ptr=&B[7*len];

    register int b0_reg, b1_reg, b2_reg, b3_reg, b4_reg, b5_reg, b6_reg, b7_reg;

    for(int i=0; i<len; i++)
    {
        a0_reg=A[i];
        a1_reg=A[len+i];
        a2_reg=A[2*len+i];
        a3_reg=A[3*len+i];

        b0_reg=*b0_ptr++;
        b1_reg=*b1_ptr++;
        b2_reg=*b2_ptr++;
        b3_reg=*b3_ptr++;
        b4_reg=*b4_ptr++;
        b5_reg=*b5_ptr++;
        b6_reg=*b6_ptr++;
        b7_reg=*b7_ptr++;

        c00_reg+=a0_reg*b0_reg;
        c01_reg+=a0_reg*b1_reg;
        c02_reg+=a0_reg*b2_reg;
        c03_reg+=a0_reg*b3_reg;
        c04_reg+=a0_reg*b4_reg;
        c05_reg+=a0_reg*b5_reg;
        c06_reg+=a0_reg*b6_reg;
        c07_reg+=a0_reg*b7_reg;

        c10_reg+=a1_reg*b0_reg;
        c11_reg+=a1_reg*b1_reg;
        c12_reg+=a1_reg*b2_reg;
        c13_reg+=a1_reg*b3_reg;
        c14_reg+=a1_reg*b4_reg;
        c15_reg+=a1_reg*b5_reg;
        c16_reg+=a1_reg*b6_reg;
        c17_reg+=a1_reg*b7_reg;

        c20_reg+=a2_reg*b0_reg;
        c21_reg+=a2_reg*b1_reg;
        c22_reg+=a2_reg*b2_reg;
        c23_reg+=a2_reg*b3_reg;
        c24_reg+=a2_reg*b4_reg;
        c25_reg+=a2_reg*b5_reg;
        c26_reg+=a2_reg*b6_reg;
        c27_reg+=a2_reg*b7_reg;

        c30_reg+=a3_reg*b0_reg;
        c31_reg+=a3_reg*b1_reg;
        c32_reg+=a3_reg*b2_reg;
        c33_reg+=a3_reg*b3_reg;
        c34_reg+=a3_reg*b4_reg;
        c35_reg+=a3_reg*b5_reg;
        c36_reg+=a3_reg*b6_reg;
        c37_reg+=a3_reg*b7_reg;
    }

    result[0]=c00_reg;
    result[1]=c01_reg;
    result[2]=c02_reg;
    result[3]=c03_reg;
    result[4]=c04_reg;
    result[5]=c05_reg;
    result[6]=c06_reg;
    result[7]=c07_reg;

    result[len+0]=c10_reg;
    result[len+1]=c11_reg;
    result[len+2]=c12_reg;
    result[len+3]=c13_reg;
    result[len+4]=c14_reg;
    result[len+5]=c15_reg;
    result[len+6]=c16_reg;
    result[len+7]=c17_reg;

    result[2*len+0]=c20_reg;
    result[2*len+1]=c21_reg;
    result[2*len+2]=c22_reg;
    result[2*len+3]=c23_reg;
    result[2*len+4]=c24_reg;
    result[2*len+5]=c25_reg;
    result[2*len+6]=c26_reg;
    result[2*len+7]=c27_reg;

    result[3*len+0]=c30_reg;
    result[3*len+1]=c31_reg;
    result[3*len+2]=c32_reg;
    result[3*len+3]=c33_reg;
    result[3*len+4]=c34_reg;
    result[3*len+5]=c35_reg;
    result[3*len+6]=c36_reg;
    result[3*len+7]=c37_reg;

    return;
}

int main()
{
    using namespace std;

    input(matA, matB);  //随机生成A和B
    int *matC = matCm, *matC2 = matCm2, n = 1;
    memcpy(matC, matA, sizeof(int[N*N]));

    cout<<"start main"<<endl;
    
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
