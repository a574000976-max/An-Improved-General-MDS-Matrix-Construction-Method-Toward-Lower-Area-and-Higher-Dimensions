#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <bitset>
#include <iostream>
using namespace std;
#ifndef SIZE
#error "Compile with -DSIZE=<binary matrix dimension>"
#endif
#define MULTI_THREAD_FLAG 0
#define THREAD_NUM 1
typedef bitset<SIZE> ROW;
struct xpair { int src=0; int dst=0; bool flag=false; };
struct thread_data { vector<xpair> seq; int gap; int start; int len; };
vector<ROW> get_matrix();
#endif
