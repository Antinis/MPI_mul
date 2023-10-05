// original author: Xu Qiyuan
#pragma once
#include <omp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>

#define __X_N 10001
int __X_matA[__X_N*__X_N], __X_matB[__X_N*__X_N],
    __X_tv_m[__X_N], __X_tv_m2[__X_N];
struct timespec __X_begin, __X_end;

//保证只构造一次
void input(int* matA, int* matB) {
    static int pid = -1;
    if (pid != -1) {
        puts("input can only be called once by one process/thread");
        abort();
    }
    pid = 0;

    srand(time(NULL));

// 生成随机数，首先生成在全局变量里，然后复制到传入指针中
#pragma omp parallel for
    for (int i=0;i<__X_N*__X_N;++i)
        __X_matA[i] = rand();
#pragma omp parallel for
    for (int i=0;i<__X_N*__X_N;++i)
        __X_matB[i] = rand();

#pragma omp parallel for
    for (int i = 0; i < __X_N ; ++i)
        for (int j = 0; j < __X_N ; ++j)
            matA[i*__X_N + j] = __X_matA[j*__X_N + i];
#pragma omp parallel for
    for (int i = 0; i < __X_N ; ++i)
        for (int j = 0; j < __X_N ; ++j)
            matB[i*__X_N + j] = __X_matB[j*__X_N + i];
    
    // std::cout<<"Timer start"<<std::endl;
    clock_gettime(CLOCK_MONOTONIC, &__X_begin);
}

void output(int* result, int n)
{
    // using namespace std;

    // std::cout<<"Timer end"<<std::endl;
    clock_gettime(CLOCK_MONOTONIC, &__X_end);
    printf("Performance : %lf Gops\n", (double)__X_N * __X_N * (__X_N + 1) * n / 1000000000/ ((double)__X_end.tv_sec - __X_begin.tv_sec + 0.000000001 * (__X_end.tv_nsec - __X_begin.tv_nsec)));
    // cout<<__X_end.tv_sec - __X_begin.tv_sec;

    srand(__X_end.tv_nsec);
    int test_on = rand() % __X_N,
        *__X_tv = __X_tv_m, *__X_tv2 = __X_tv_m2;
    // test_on=1500;
    for (int i=0;i<__X_N;++i)
        __X_tv[i] = __X_matA[i*__X_N + test_on];    //随机选mat_A的一列

    // cout<<"check calc start"<<endl;
    // cout<<"check row: "<<test_on<<endl;
    for (int k=0;k<n;++k) //阶数（默认1阶，n=1）
    {
        #pragma omp parallel for
        for (int i=0;i<__X_N*__X_N;++i)
            __X_matA[i] += __X_matB[i];

        #pragma omp parallel for
        for (int i=0;i<__X_N;++i)   //行号
        {
            int sum = 0, row = i*__X_N;
            for (int k=0;k<__X_N;++k) 
            {
                sum += __X_tv[k] * __X_matA[row + k];
            }
            __X_tv2[i] = sum;
        }
        int *t = __X_tv; __X_tv = __X_tv2; __X_tv2 = t;
    }
    // cout<<"check calc over"<<endl;

    //随机取答案的一列验证正确性
    for (int i=0;i<__X_N;++i)
        if (__X_tv[i] != result[test_on * __X_N + i])
        {
            printf("Result is incorrect : should be %d but actually %d at %d,%d\n", __X_tv[i], result[test_on * __X_N + i], test_on, i);
            abort();
        }
    
    puts("Validation passed.");
}
