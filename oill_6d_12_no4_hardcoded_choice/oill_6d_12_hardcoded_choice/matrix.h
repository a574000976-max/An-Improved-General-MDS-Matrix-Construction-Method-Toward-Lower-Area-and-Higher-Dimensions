#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <bitset>
#include <iostream>

using namespace std;

#define MULTI_THREAD_FLAG 1
#define THREAD_NUM 4

// 6-dimensional L-power matrices only.
// 4-bit entries have been removed.
// Change CHOICE from 1 to 12, then run: make clean && make && ./evaluate
#define CHOICE 12

#if CHOICE == 1
    #define SIZE 8
    #define FILENAME "6x6_n8_L_inv_pow_2_result.txt"

#elif CHOICE == 2
    #define SIZE 8
    #define FILENAME "6x6_n8_L_inv_pow_1_result.txt"

#elif CHOICE == 3
    #define SIZE 8
    #define FILENAME "6x6_n8_L_pow_1_result.txt"

#elif CHOICE == 4
    #define SIZE 16
    #define FILENAME "6x6_n16_L_inv_pow_2_result.txt"

#elif CHOICE == 5
    #define SIZE 16
    #define FILENAME "6x6_n16_L_inv_pow_1_result.txt"

#elif CHOICE == 6
    #define SIZE 16
    #define FILENAME "6x6_n16_L_pow_1_result.txt"

#elif CHOICE == 7
    #define SIZE 32
    #define FILENAME "6x6_n32_L_inv_pow_2_result.txt"

#elif CHOICE == 8
    #define SIZE 32
    #define FILENAME "6x6_n32_L_inv_pow_1_result.txt"

#elif CHOICE == 9
    #define SIZE 32
    #define FILENAME "6x6_n32_L_pow_1_result.txt"

#elif CHOICE == 10
    #define SIZE 64
    #define FILENAME "6x6_n64_L_inv_pow_2_result.txt"

#elif CHOICE == 11
    #define SIZE 64
    #define FILENAME "6x6_n64_L_inv_pow_1_result.txt"

#elif CHOICE == 12
    #define SIZE 64
    #define FILENAME "6x6_n64_L_pow_1_result.txt"

#else
    #error "Unsupported CHOICE. Please set CHOICE from 1 to 12."
#endif

typedef bitset<SIZE> ROW;

typedef struct{
    int src;
    int dst;
    bool flag;
}xpair;

typedef struct
{
    vector<xpair> seq;
    int gap;
    int start;
    int len;
}thread_data;

vector<ROW> get_matrix();

#endif
