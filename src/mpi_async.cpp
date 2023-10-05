#include <iostream>
#include "hw.h"
#include <omp.h>
#include <mpi.h>
#include <cstdlib>
#include <time.h>
#include <string.h>

#define N __X_N //10001

int matA[N*N], matB[N*N],
    matCm[N*N], matCm2[N*N];

int RecvTmp[3][N*N/3];
clock_t send_begin, send_end;
clock_t recv_begin, recv_end;

void calc_4x4(int *result, int *A, int *B, int slave_len, int N_len)
{
    //B转置过
    register int c00_reg=0, c01_reg=0, c02_reg=0, c03_reg=0;
    register int c10_reg=0, c11_reg=0, c12_reg=0, c13_reg=0;
    register int c20_reg=0, c21_reg=0, c22_reg=0, c23_reg=0;
    register int c30_reg=0, c31_reg=0, c32_reg=0, c33_reg=0;

    int *b0_ptr=&B[0];
    int *b1_ptr=&B[N_len];
    int *b2_ptr=&B[2*N_len];
    int *b3_ptr=&B[3*N_len];

    int *a0_ptr=&A[0];
    int *a1_ptr=&A[N_len];
    int *a2_ptr=&A[2*N_len];
    int *a3_ptr=&A[3*N_len];

    register int a0_reg, a1_reg, a2_reg, a3_reg;
    register int b0_reg, b1_reg, b2_reg, b3_reg;

    for(int i=0; i<N_len; i++)
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

    result[slave_len+0]=c10_reg;
    result[slave_len+1]=c11_reg;
    result[slave_len+2]=c12_reg;
    result[slave_len+3]=c13_reg;

    result[2*slave_len+0]=c20_reg;
    result[2*slave_len+1]=c21_reg;
    result[2*slave_len+2]=c22_reg;
    result[2*slave_len+3]=c23_reg;

    result[3*slave_len+0]=c30_reg;
    result[3*slave_len+1]=c31_reg;
    result[3*slave_len+2]=c32_reg;
    result[3*slave_len+3]=c33_reg;

    return;
}

void calc_4x4_master(int *result, int *A, int *B, int len)
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

void calc_1x4(int *result, int *A, int *B, int len)
{
    //B转置过
    register int a_reg;
    register int c0_reg, c1_reg, c2_reg, c3_reg;
    c0_reg=0;
    c1_reg=0;
    c2_reg=0;
    c3_reg=0;

    int *b0_ptr=&B[0];
    int *b1_ptr=&B[len];
    int *b2_ptr=&B[2*len];
    int *b3_ptr=&B[3*len];

    for(int i=0; i<len; i++)
    {
        a_reg=A[i];

        c0_reg+=a_reg* *b0_ptr++;
        c1_reg+=a_reg* *b1_ptr++;
        c2_reg+=a_reg* *b2_ptr++;
        c3_reg+=a_reg* *b3_ptr++;
    }

    result[0]=c0_reg;
    result[1]=c1_reg;
    result[2]=c2_reg;
    result[3]=c3_reg;

    return;
}

int kernel_1x1(int *A, int *B, int len)
{
    register int ans=0;
    int *a_ptr=A;
    int *b_ptr=B;
    for(int i=0; i<len; i++)
        ans+=(*a_ptr++)*(*b_ptr++);
    
    return ans;
}

int main()
{
    using namespace std;

    int n = 8;
    int slave_len=(N/8)*4; // 4的倍数


    MPI_Init(NULL, NULL);

    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if(rank==0) // master管理节点，同时作为slave 1计算节点
    {
        cout<<"master: preparing..."<<endl;

        input(matA, matB);  //随机生成A和B

        int *matC = matCm, *matC2 = matCm2;
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

        // MPI_Send(数据指针, 数据长度, 数据类型, 目标节点, TAG, MPI_COMM_WORLD);
        // MPI_Recv(数据指针, 数据长度, 数据类型, 来源节点, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // master节点自适应任意形状矩阵
        for (int k=0;k<n;++k)
        {
            cout<<"master: processing "<<"n: "<<n<<"; k: "<<k<<endl;

            #pragma omp parallel for
            for (int i=0;i<N*N;++i)
                matA[i] += matB[i];

            // 下发数据到节点，有四个slave节点
            cout<<"master: sending..."<<endl;
            send_begin=clock();
            MPI_Request request1;
            for(int slave=1; slave<=world_size-1; slave++)
            {
                // 左矩阵发送
                MPI_Isend(&matC[(slave)/2*N*slave_len], N*slave_len, MPI_INT, slave, 0, MPI_COMM_WORLD, &request1);
                // 右矩阵发送
                MPI_Isend(&matA[(slave)%2*N*slave_len], N*slave_len, MPI_INT, slave, 1, MPI_COMM_WORLD, &request1);
            }
            send_end=clock();
            cout<<"master: sending over"<<endl;
            cout<<"master: sending time usage: "<<(double)(send_end-send_begin)/CLOCKS_PER_SEC*1000<<" ms"<<endl;

            //master节点开始自身的计算，作为slave 0
            cout<<"slave "<<rank<<": calc start..."<<endl;
            #pragma omp parallel for
            for (int i=0; i<slave_len; i+=4)
                for (int j=0; j<slave_len; j+=4)
                    calc_4x4_master(&matC2[i*N+j], &matC[i*N], &matA[j*N], N);
            
            cout<<"slave "<<rank<<": calc over"<<endl;
            
            // 同步发送
            MPI_Wait(&request1, MPI_STATUSES_IGNORE);

            //接受来自计算节点的结果
            cout<<"master: recieving..."<<endl;
            recv_begin=clock();
            MPI_Request request2[3];
            MPI_Status status[3];
            MPI_Irecv(RecvTmp[1], slave_len*slave_len, MPI_INT, 2, 2, MPI_COMM_WORLD, &request2[0]);
            MPI_Irecv(RecvTmp[0], slave_len*slave_len, MPI_INT, 1, 2, MPI_COMM_WORLD, &request2[1]);
            MPI_Irecv(RecvTmp[2], slave_len*slave_len, MPI_INT, 3, 2, MPI_COMM_WORLD, &request2[2]);

            recv_end=clock();
            cout<<"master: recieving over"<<endl;
            cout<<"master: sending time usage: "<<(double)(recv_end-recv_begin)/CLOCKS_PER_SEC*1000<<" ms"<<endl;

            // 计算边角，下方
            for(int i=slave_len*2; i<N; i++)
                #pragma omp parallel for
                for(int j=0; j<slave_len*2; j+=4)
                    calc_1x4(&matC2[i*N+j], &matC[i*N], &matA[j*N], N);
            
            // 计算边角，右方
            #pragma omp parallel for
            for(int i=0; i<N; i++)
                for(int j=slave_len*2; j<N; j++)
                    matC2[i*N+j]=kernel_1x1(&matC[i*N], &matA[j*N], N);
            
            // 同步接收
            MPI_Waitall(3, request2, MPI_STATUS_IGNORE);

            for(int slave=1; slave<world_size; slave++)
            {
                int bias_i=(slave)/2*slave_len;
                int bias_j=(slave)%2*slave_len;

                #pragma omp parallel for
                for(int i=0; i<slave_len; i++)
                    for(int j=0; j<slave_len; j++)
                        matC2[(i+bias_i)*N+j+bias_j]=RecvTmp[slave-1][i*slave_len+j];
            }
            
            int *t = matC; matC = matC2; matC2 = t;
        }

        cout<<"master: calc over"<<endl;
        output(matC, n);
    }
    else    // slave节点
    {
        // slave不需要知道自己的绝对位置
        cout<<"slave: rank: "<<rank<<endl;

        int *matC = matCm, *matC2 = matCm2;

        // MPI_Send(数据指针, 数据长度, 数据类型, 目标节点, TAG, MPI_COMM_WORLD);
        // MPI_Recv(数据指针, 数据长度, 数据类型, 来源节点, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int k=0;k<n;++k)
        {
            // 左矩阵接收
            cout<<"slave "<<rank<<": recieving..."<<endl;
            MPI_Recv(matC, N*slave_len, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // 右矩阵接收，转置
            MPI_Recv(matA, N*slave_len, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            cout<<"slave "<<rank<<": recieving over"<<endl;
            
            // 计算
            cout<<"slave "<<rank<<": calc start..."<<endl;
            #pragma omp parallel for
            for (int i=0; i<slave_len; i+=4)
                for (int j=0; j<slave_len; j+=4)
                {
                    //计算第i行第j列的结果
                    calc_4x4(&matC2[i*slave_len+j], &matC[i*N], &matA[j*N], slave_len, N);
                }
            cout<<"slave "<<rank<<": calc over"<<endl;
            
            // 计算结果传回master节点
            cout<<"slave "<<rank<<": sending matC2..."<<endl;
            MPI_Send(matC2, slave_len*slave_len, MPI_INT, 0, 2, MPI_COMM_WORLD);
            cout<<"slave "<<rank<<": sending matC2 over"<<endl;

            // cout<<matC2[slave_len*100]<<endl;
        }
    }

    MPI_Finalize();

    return 0;
}
