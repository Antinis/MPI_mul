#include <iostream>
#include "hw.h"
#include <omp.h>
#include <string.h>
#include <immintrin.h>

#pragma pack(4)

#define N __X_N //10001

int matA[N*N], matB[N*N],
    matCm[N*N], matCm2[N*N];

void calc_1x8(int *result, int *A, int *B, int len)
{
    //A没有转置，B没有转置

    __m256i vec_a0, vec_a1, vec_a2, vec_a3;
    __m256i vec_c0, vec_c1, vec_c2, vec_c3;
    __m256i vec_b;

    __m256i vec_tmp0, vec_tmp1, vec_tmp2, vec_tmp3;

    vec_c0=_mm256_set1_epi32(0);
    vec_c1=_mm256_set1_epi32(0);
    vec_c2=_mm256_set1_epi32(0);
    vec_c3=_mm256_set1_epi32(0);

    int *a0_ptr=&A[0];
    int *a1_ptr=&A[len];
    int *a2_ptr=&A[2*len];
    int *a3_ptr=&A[3*len];

    for(int i=0; i<len; i++)
    {
        vec_b=_mm256_load_si256(reinterpret_cast<__m256i const*>(B+i*len));

        vec_a0=_mm256_set1_epi32(*a0_ptr++);
        vec_a1=_mm256_set1_epi32(*a1_ptr++);
        vec_a2=_mm256_set1_epi32(*a2_ptr++);
        vec_a3=_mm256_set1_epi32(*a3_ptr++);
        
        vec_tmp0=_mm256_mullo_epi32(vec_a0, vec_b);
        vec_tmp1=_mm256_mullo_epi32(vec_a1, vec_b);
        vec_tmp2=_mm256_mullo_epi32(vec_a2, vec_b);
        vec_tmp3=_mm256_mullo_epi32(vec_a3, vec_b);

        vec_c0=_mm256_add_epi32(vec_c0, vec_tmp0);
        vec_c1=_mm256_add_epi32(vec_c1, vec_tmp0);
        vec_c2=_mm256_add_epi32(vec_c2, vec_tmp0);
        vec_c3=_mm256_add_epi32(vec_c3, vec_tmp0);
    }

    _mm256_store_si256(reinterpret_cast<__m256i *>(result), vec_c0);
    _mm256_store_si256(reinterpret_cast<__m256i *>(result+len), vec_c1);
    _mm256_store_si256(reinterpret_cast<__m256i *>(result+2*len), vec_c2);
    _mm256_store_si256(reinterpret_cast<__m256i *>(result+3*len), vec_c3);

    return;
}

int main()
{
    using namespace std;

    input(matA, matB);  //随机生成A和B
    int *matC = matCm, *matC2 = matCm2, n = 1;
    memcpy(matC, matA, sizeof(int[N*N]));

    // //转置A
    // #pragma omp parallel for schedule(guided, 100)
    // for (int i=0;i<N;++i)
    //     for (int j=0;j<i;++j) {
    //         int t = matA[i*N+j];
    //         matA[i*N+j] = matA[j*N+i];
    //         matA[j*N+i] = t;
    //     }

    // //转置B
    // #pragma omp parallel for schedule(guided, 100)
    // for (int i=0;i<N;++i)
    //     for (int j=0;j<i;++j) {
    //         int t = matB[i*N+j];
    //         matB[i*N+j] = matB[j*N+i];
    //         matB[j*N+i] = t;
        // }

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
                //matC没有转置，matA没有转置
                calc_1x8(&matC2[i*N+j], &matC[i*N], &matA[j], N);
            }

        int *t = matC; matC = matC2; matC2 = t;
    }

    cout<<"calc over"<<endl;
    output(matC, n);
}
