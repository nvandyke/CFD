
#include "Structures.h"

#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstring>

using namespace std;

vec::vec() {
	sizeI = 1;
	sizeJ = 1;
	totalSize = 1;
	data = new double[totalSize];
	setZero();

}


vec::vec(int i, int j) {
	sizeI = i;
	sizeJ = j;
	totalSize = sizeI * sizeJ;
	data = new double[totalSize];
	setZero();
}

vec::vec(const vec& rhs) {
	
	sizeI = rhs.rows();
	sizeJ = rhs.cols();
	totalSize = sizeI * sizeJ;
	if (data) {
		delete[] data;
	}
	data = new double[totalSize];
	std::memcpy(data, rhs.data, rhs.totalSize * sizeof(double));
	/*for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ;j++){
			setAt(i, j, v.getAt(i, j));
		}
	}*/
}

vec::~vec() {
	
}

double vec::dot(const vec& v) {
	//ensure the same size
	assert(sizeCheck(*this, v));
	assert(rows() == 1 || cols() == 1);
	double value = 0.0;
	for (int i = 0; i < rows(); i++) {
		for (int j = 0; j < cols(); j++) {
			value += this->getAt(i, j) * v.getAt(i, j);
		}
	}

	return value;
}

void vec::print(ostream& os){
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			os << this->getAt(i,j) << " ";
		}
		os << endl;
	}
}

void vec::setAt(int i, int j, double df){
	if ((i >= 0) && (j >= 0) && (i <= sizeI) && (j <= sizeJ)) {
		int index = i * int(sizeJ) + j;
		data[index] = df;
	} else {
		cout << "out of bounds" << endl;
		assert(false);
	}
}

const double vec::getAt(int i, int j) const {
	if ((i >= 0) && (j >= 0) && (i <= sizeI) && (j <= sizeJ)) {
		int index = i * int(sizeJ) + j;
		return data[index];
	} else {
		assert(false);
		return 0;
	}
}

vec& vec::operator+(const vec& b){
	assert(sizeCheck(*this, b));
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) + b.getAt(i, j));

		}
	}
	return *this;
}

vec& vec::operator-(const vec& b){
	assert(sizeCheck(*this, b));
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) - b.getAt(i, j));

		}
	}
	return *this;
}

vec& vec::operator*(const vec& b){
	assert(sizeCheck(*this, b));
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) * b.getAt(i, j));

		}
	}
	return *this;
}

vec& vec::operator/(const vec& b){
	assert(sizeCheck(*this, b));
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) / b.getAt(i, j));

		}
	}
	return *this;
}

vec& vec::operator=(const vec& rhs) {
	//assert(sizeCheck(*this, b));
	if (&rhs == this) {
		return *this;
	} 

	if (totalSize == rhs.totalSize) {
		std::memcpy(data, rhs.data, rhs.totalSize * sizeof(double));
	} else {
		delete[] data;
		data = new double[rhs.totalSize]();
		std::memcpy(data, rhs.data, rhs.totalSize * sizeof(double));
	}
	sizeI = rhs.sizeI;
	sizeJ = rhs.sizeJ;
	totalSize = rhs.totalSize;

	/*
	double* data2 = new double[sizeI * sizeJ];

	std::copy(data2, data2 + sizeI * sizeJ, b.data);

	delete[] data;
	data = data2;
	*/
	/*for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			v.setAt(i, j, b.getAt(i, j));

		}
	}*/
	return *this;
}


vec& vec::operator+(const double& b) {
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) + b);

		}
	}
	return *this;
}

vec& vec::operator-(const double& b) {
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) - b);

		}
	}
	return *this;
}

vec& vec::operator*(const double& b) {
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) * b);

		}
	}
	return *this;
}

vec& vec::operator/(const double& b) {
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			this->setAt(i, j, this->getAt(i, j) / b);

		}
	}
	return *this;
}

vec vec::transpose() {
	vec v(sizeJ, sizeI);
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			v.setAt(j, i, this->getAt(i, j));

		}
	}
	return v;
}

double vec::max() {
	double max = getAt(0, 0);
	double value;
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			value = this->getAt(i, j);
			if (value > max)
				max = value;

		}
	}
	return max;
}

double vec::min() {
	double min = getAt(0, 0);
	double value;
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			value = this->getAt(i, j);
			if (value < min)
				min = value;

		}
	}
	return min;
}

const vec vec::getBlock(int i, int j, int x, int y) const {
	vec v(x, y);
	int c = 0;
	for (int a = i; a < i + x; a++,c++) {
		int d = 0;
		for (int b = j; b < j + y; b++,d++) {
			v.setAt(c, d, this->getAt(a, b));
		}
	}
	return v;
}

void vec::setBlock(int i, int j, vec v) {
	int x = int(v.sizeI);
	int y = int(v.sizeJ);
	int c = 0;
	for (int a = i; a < i + x; a++,c++) {
		int d = 0;
		for (int b = j; b < j + y; b++,d++) {
			this->setAt(a, b, v.getAt(c, d));
		}
	}
}

int vec::rows() const {
	return int(sizeI);
}

int vec::cols() const {
	return int(sizeJ);
}

vec sqrt(vec v) {
	vec out(v.rows(), v.cols());
	for (int i = 0; i < v.rows(); i++) {
		for (int j = 0; j < v.cols(); j++) {
			out.setAt(i, j, sqrt(v.getAt(i, j)));
		}
	}
	return out;
}
double vec::norm() {
	double value = 0;
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			value += this->getAt(i, j) * this->getAt(i, j);
		}
	}
	return sqrt(value);
}

void vec::resize(int i, int j) {
	delete[] data;
	sizeI = i;
	sizeJ = j;
	totalSize = sizeI * sizeJ;
	data = new double[totalSize];
	setZero();
}


bool sizeCheck(const vec& v1, const vec& v2) {
	int i1 = v1.rows();
	int i2 = v2.rows();
	int j1 = v1.cols();
	int j2 = v2.cols();

	if (i1 == i2) {
		if (j1 == j2) {
			return true;
		}
	}
	return false;
}

vec vec::matMult(vec b) {
	int n = this->rows();
	int m = this->cols();
	int a = b.rows();
	int p = b.cols();

	assert(m == a);

	vec v(n,p);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p; j++) {
			double sum = 0;
			for (int k = 0; k < m; k++) {
				sum += this->getAt(i, k) * b.getAt(k, j);
			}
			v.setAt(i, j, sum);
		}
	}
	return v;

}

void vec::setZero() {
	for (int i = 0; i < this->rows(); i++) {
		for (int j = 0; j < this->cols(); j++) {
			this->setAt(i, j, 0);
		}
	}

}
