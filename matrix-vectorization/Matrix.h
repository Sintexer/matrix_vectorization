#pragma once

#define MT_M 8
#define MT_N 4
#define MT_K 5
#define CELL_ROWS 4
#define CELL_COLS 4

class Matrix {
protected:
    float **data{};
    int rows{};
    int cols{};

    const float MIN_VALUE = 0.;
    const float MAX_VALUE = 100.;

public:
    Matrix(int rows_, int cols_) : rows(rows_), cols(cols_) {
        data = allocate(rows, cols);
        fill();
    }

    Matrix() : rows(MT_M), cols(MT_N) {
        data = allocate(MT_M, MT_N);
        fill();
    }

    Matrix(const Matrix &m);

    ~Matrix();

    static Matrix createInnerResultMatrix() {
        return Matrix(MT_N, MT_K);
    }

    static Matrix multiply(Matrix& m1, Matrix& m2);

    Matrix operator+(Matrix &matrix1);

    Matrix operator*(Matrix &matrix1);

    Matrix& operator=(const Matrix &matrix);

    void add(Matrix& m);

    const void print();

    bool operator==(const Matrix &rhs) const {
        return data == rhs.data &&
               rows == rhs.rows &&
               cols == rhs.cols;
    }

    bool operator!=(const Matrix &rhs) const {
        return !(rhs == *this);
    }

    float **getData() const {
        return data;
    }

    int getRows() const {
        return rows;
    }

    int getCols() const {
        return cols;
    }

    void fill();

    void fill(float value);

protected:
    void clear();

    static float **allocate(int rows, int columns);

    float getRandomFloat(float min, float max);

    void fill(float **matrix);

};

