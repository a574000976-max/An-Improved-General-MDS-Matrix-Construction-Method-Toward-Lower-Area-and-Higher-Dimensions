#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <bitset>
#include <iostream>

using namespace std;

#define MULTI_THREAD_FLAG 1
#define THREAD_NUM 4

// Only 15 choices are kept for the 9-dimensional L-power matrices.
// Change this value from 1 to 15, then run: make clean && make && ./evaluate
#define CHOICE 15

#if CHOICE == 1
    #define SIZE 16
    #define FILENAME "GF2_16_L_inv_pow_2.txt"

#elif CHOICE == 2
    #define SIZE 16
    #define FILENAME "GF2_16_L_inv_pow_1.txt"

#elif CHOICE == 3
    #define SIZE 16
    #define FILENAME "GF2_16_L_pow_1.txt"

#elif CHOICE == 4
    #define SIZE 16
    #define FILENAME "GF2_16_L_pow_2.txt"

#elif CHOICE == 5
    #define SIZE 16
    #define FILENAME "GF2_16_L_pow_3.txt"

#elif CHOICE == 6
    #define SIZE 32
    #define FILENAME "GF2_32_L_inv_pow_2.txt"

#elif CHOICE == 7
    #define SIZE 32
    #define FILENAME "GF2_32_L_inv_pow_1.txt"

#elif CHOICE == 8
    #define SIZE 32
    #define FILENAME "GF2_32_L_pow_1.txt"

#elif CHOICE == 9
    #define SIZE 32
    #define FILENAME "GF2_32_L_pow_2.txt"

#elif CHOICE == 10
    #define SIZE 32
    #define FILENAME "GF2_32_L_pow_3.txt"

#elif CHOICE == 11
    #define SIZE 64
    #define FILENAME "GF2_64_L_inv_pow_2.txt"

#elif CHOICE == 12
    #define SIZE 64
    #define FILENAME "GF2_64_L_inv_pow_1.txt"

#elif CHOICE == 13
    #define SIZE 64
    #define FILENAME "GF2_64_L_pow_1.txt"

#elif CHOICE == 14
    #define SIZE 64
    #define FILENAME "GF2_64_L_pow_2.txt"

#elif CHOICE == 15
    #define SIZE 64
    #define FILENAME "GF2_64_L_pow_3.txt"

#else
    #error "Unsupported CHOICE. Please set CHOICE from 1 to 15."
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
