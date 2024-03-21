/*
 * matrix.cpp
 */

#include <stdexcept>
#include "block.h"

#define EPS 1e-10

using std::ostream;  using std::istream;  using std::endl;
using std::domain_error;

/* PUBLIC MEMBER FUNCTIONS
 ********************************/

Block::Block(int rows, int cols): rows_(rows), cols_(cols) {
    allocSpace();
    for (int i = 0; i < size(); ++i) {
        p[i] = 0;
    }

}

Block::Block(int* a, int rows, int cols): rows_(rows), cols_(cols) {
    allocSpace();
    for (int i = 0; i < size(); ++i) {
        p[i] = a[i];
    }
}

Block::Block(): rows_(1), cols_(1) {
    allocSpace();
    p[0] = 0;
}

Block::~Block() {
    delete[] p;
}

Block::Block(const Block& m): rows_(m.rows_), cols_(m.cols_) {
    allocSpace();
    for (int i = 0; i < size(); ++i) {
        p[i] = m.p[i];
    }
}

Block& Block::operator=(const Block& m) {
    if (this == &m) {
        return *this;
    }

    if (rows_ != m.rows_ || cols_ != m.cols_) {
        delete[] p;

        rows_ = m.rows_;
        cols_ = m.cols_;
        allocSpace();
    }

    for (int i = 0; i < size(); ++i) {
        p[i] = m.p[i];
    }
    return *this;
}

Block& Block::operator+=(const Block& m) {
    for (int i = 0; i < size(); ++i) {
        p[i] += m.p[i];
    }
    return *this;
}

Block& Block::operator-=(const Block& m) {
    for (int i = 0; i < size(); ++i) {
        p[i] -= m.p[i];
    }
    return *this;
}
/*
Block& Block::operator*=(const Block& m) {
    Block temp(rows_, m.cols_);
    Block mT = m.transpose();
    for (int i = 0; i < temp.rows_; ++i) {
        for (int j = 0; j < temp.cols_; ++j) {
            for (int k = 0; k < cols_; ++k) {
                temp.p[i*cols_+j] += (p[i*cols_+k] * m.p[j*m.cols_+k]);
            }
        }
    }
    return (*this = temp);
}
*/
Block& Block::operator*=(int num) {
    for (int i = 0; i < size(); ++i) {
        p[i] *= num;
    }
    return *this;
}

Block& Block::operator/=(int num) {
    for (int i = 0; i < size(); ++i) {
        p[i] /= num;
    }
    return *this;
}


Block Block::transpose() {
    Block ret(cols_, rows_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            ret.p[j * cols_ + i] = p[i * cols_ + j];
        }
    }
    return ret;
}

/* PRIVATE HELPER FUNCTIONS
 ********************************/

void Block::allocSpace() {
    p = new int [size()];
}



/* NON-MEMBER FUNCTIONS
 ********************************/

Block operator+(const Block& m1, const Block& m2) {
    Block temp(m1);
    return (temp += m2);
}

Block operator-(const Block& m1, const Block& m2) {
    Block temp(m1);
    return (temp -= m2);
}
/*
Block operator*(const Block& m1, const Block& m2) {
    Block temp(m1);
    return (temp *= m2);
}
*/
Block operator*(const Block& m, int num) {
    Block temp(m);
    return (temp *= num);
}

Block operator*(int num, const Block& m) {
    return (m * num);
}

Block operator/(const Block& m, int num) {
    Block temp(m);
    return (temp /= num);
}

ostream& operator<<(ostream& os, const Block& m) {
    for (int i = 0; i < m.rows_; ++i) {
        os << m.p[i];
        for (int j = 1; j < m.cols_; ++j) {
            os << " " << m.p[i * m.cols_ + j];
        }
        os << endl;
    }
    return os;
}

istream& operator>>(istream& is, Block& m) {
    for (int i = 0; i < m.size(); ++i) {
        is >> m.p[i];
    }
    return is;
}

int Block::rows() const {
    return rows_;
}

int Block::cols() const {
    return cols_;
}

void Block::print(std::ostream& os) {
    int index = 0;
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            os << p[i * cols_ + j] << " ";
        }
        os << endl;
    }
}

void Block::setBlock(int i, int j, Block& m) {
    for (int a = i, c = 0; a < i + m.rows(); a++, c++) {
        for (int b = j, d = 0; b < j + m.cols(); b++, d++) {
            p[a * cols_ + b] = m.p[c * m.cols_ + d];
        }
    }

}

const Block Block::getBlock(int i, int j, int x, int y) const {
    Block m(x, y);
    for (int a = i, c = 0; a < i + x; a++, c++) {
        for (int b = j, d = 0; b < j + y; b++, d++) {
            m.p[c * m.cols_ + d] = p[a * cols_ + b];
        }
    }
    return m;
}

int Block::max() {
    int max = p[0];
    for (int i = 0; i < size(); i++) {
        if (p[i] > max) {
            max = p[i];
        }
    }
    return max;
}

int Block::size() {
    return rows_ * cols_;
}

