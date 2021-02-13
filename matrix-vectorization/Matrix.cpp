
#include <cstdlib>
#include <random>
#include <iostream>
#include <iomanip>
#include "Matrix.h"

Matrix::Matrix(const Matrix& matrix) {
    data = new float* [matrix.rows];
    cols = matrix.cols;
    rows = matrix.rows;
    for (int i = 0; i < MT_ROWS; ++i) {
        data[i] = new float[matrix.cols];
        for (int j = 0; j < MT_COLS; ++j) {
            data[i][j] = matrix.getData()[i][j];
        }
    }
}

Matrix::~Matrix() {
    Matrix::clear();
}

Matrix Matrix::multiply(Matrix& m1, Matrix& m2)
{
    auto data = m1.data;
    auto data_ = m2.getData();
    auto cols = m1.cols;
    auto cols_ = m2.getCols();
    auto rows = m1.rows;
    auto resultData = allocate(rows, cols_);

    for (int i = 0; i < MT_ROWS; ++i) {
        auto resultRow = resultData[i];
        auto row = data[i];
        for (int j = 0; j < MT_COLS; ++j) {
            auto val = row[j];
            auto col = data_[j];
            for (int k = 0; k < cols_; ++k) {
                resultRow[k] += val * col[k];
            }
        }
    }

    Matrix result(rows, cols_);
    result.fill(resultData);

    return result;
}

float **Matrix::allocate(int rows, int cols) {
    auto data_ = new float *[rows];
    for (int i = 0; i < rows; ++i) {
        data_[i] = new float[cols];
        for (int j = 0; j < cols; ++j)
            data_[i][j] = 0;
    }
    return data_;
}

void Matrix::fill() {
    for (int i = 0; i < MT_ROWS; ++i) {
        for (int j = 0; j < MT_COLS; ++j) {
            data[i][j] = getRandomFloat(MIN_VALUE, MAX_VALUE);
        }
    }
}

float Matrix::getRandomFloat(float min, float max) {
    float a = min + (float) (rand()) / ((float) (RAND_MAX / (max - min)));
    return a;
}

void Matrix::fill(float value) {
    for (int i = 0; i < MT_ROWS ; ++i) {
        for (int j = 0; j < MT_COLS; ++j) {
            data[i][j] = value;
        }
    }
}

void Matrix::clear() {
    for (int i = 0; i < MT_ROWS; ++i) {
        delete[] data[i];
    }
    delete data;
}

void Matrix::print() {
    for (int i = 0; i < MT_ROWS; i++) {
        for (int j = 0; j < MT_COLS; j++) {
            std::cout << std::setprecision(8) << std::setw(7) <<data[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

void Matrix::fill(float **matrix) {
    for (int i = 0; i < MT_ROWS; i++) {
        for (int j = 0; j < MT_COLS; j++) {
            data[i][j] = matrix[i][j];
        }
    }
}

Matrix Matrix::operator*(Matrix &matrix) {
    //if (cols != matrix.getRows())
    //    throw std::runtime_error("Can't multiply matrix");

    auto data_ = matrix.getData();

    auto cols_ = matrix.getCols();
    auto resultData = allocate(rows, cols_);

    for (int i = 0; i < MT_ROWS; ++i) {
        auto resultRow = resultData[i];
        auto row = data[i];
        for(int j = 0; j < MT_COLS; ++j){
            auto val = row[j];
            auto col = data_[j];
            for(int k = 0; k < cols_; ++k){
                resultRow[k] += val * col[k];
            }
        }
    }

    Matrix result(rows, cols_);
    result.fill(resultData);

    return result;

}

Matrix Matrix::operator+(Matrix &matrix1) {
    auto resultData = allocate(this->rows, this->cols);
    auto data1 = matrix1.getData();
    auto data2=data;
    for (int i=0; i<MT_ROWS; ++i) {
        auto row1=data1[i];
        auto row2=data2[i];
        for (int j=0; j< MT_COLS; ++j) {
            resultData[i][j]=row1[j]+row2[j];
        }
    }

    Matrix result(rows, cols);
    result.fill(resultData);

    return result;
}

Matrix& Matrix::operator=(const Matrix &matrix) {
    this->cols=matrix.cols;
    this->rows=matrix.rows;
    this->data=allocate(matrix.rows, matrix.cols);
    for (int i = 0; i < MT_ROWS; ++i) {
        for (int j = 0; j < MT_COLS; ++j) {
            this->data[i][j] = matrix.getData()[i][j];
        }
    }
    return *this;
}

void Matrix::add(Matrix& m)
{
    auto data1 = m.getData();
    for (int i = 0; i < MT_ROWS; ++i) {
        auto row1 = data1[i];
        for (int j = 0; j < MT_COLS; ++j) {
            data[i][j] += row1[j];
        }
    }
}

void Matrix::mult(Matrix& m)
{
    if (cols != m.getRows())
        throw std::runtime_error("Can't multiply matrix");

    auto data_ = m.getData();

    auto cols_ = m.getCols();

    for (int i = 0; i < MT_ROWS; ++i) {
        auto resultRow = data[i];
        auto row = data[i];
        for (int j = 0; j < MT_COLS; ++j) {
            auto val = row[j];
            auto col = data_[j];
            for (int k = 0; k < cols_; ++k) {
                resultRow[k] += val * col[k];
            }
        }
    }

}





