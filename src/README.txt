hw.h中为了方便调试和对比，加入了若干输出语句，并将ops单位调整为了Gops。计时点没有改变。
hw_ori.h为原始文件。

运行方法：
首先 source /opt/intel/oneapi/setvars.sh 加载MPI环境
然后salloc -N 4申请集群使用权限
然后make run运行程序
在输出中会有以Gops（使用hw.h）或者以Mops（使用hw_ori.h）为单位的算力指标。

测试Ops值：395.9257 Gops （395925.7Mops）