#include <iostream>
#include "Matrix.h"

int main() {

    const int M = 5;
    const int N = 4;
    const int K = 7;

    auto **pMatrix1 = new Matrix *[M];
    for (int i = 0; i < M; ++i) {
        pMatrix1[i] = new Matrix[N];
        for (int j = 0; j < N; ++j) {
            Matrix matrix(MT_M, MT_N);
            matrix.fill(2);
            pMatrix1[i][j]=matrix;
            //matrix.print();
        }
    }

    
    auto **pMatrix2 = new Matrix *[N];
    for (int i = 0; i < N; ++i) {
        pMatrix2[i] = new Matrix[K];
        for (int j = 0; j < K; ++j) {
            Matrix matrix(MT_N, MT_K);
            matrix.fill(2);

            pMatrix2[i][j]=matrix;
            //matrix.print();
        }
    }

    auto **result = new Matrix *[M];
    for (int i = 0; i < M; ++i) {
        result[i] = new Matrix[K];
        for (int j = 0; j < K; ++j) {
            Matrix matrix(MT_M, MT_K);
            matrix.fill(0.);
            result[i][j]=matrix;
        }
    }

    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < K; ++k) {
                Matrix m  = pMatrix1[i][j] * pMatrix2[j][k];
                result[i][k].add(m);
                if (k==K-1) {
                    result[i][k].print();
                }
            }
        }
    }



}
