#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>
using namespace std;

class vec
{
private:
	int sizeI;
	int sizeJ;
	int totalSize;
	double* data;
public:
	vec();
	vec(int i, int j);
	vec(const vec& v);
	~vec();
	double dot(const vec& v);
	void print(ostream& os = std::cout);
	void setAt(int i, int j, double df);
	const double getAt(int i, int j) const;

	//Elementwise operations vec-vec
	vec& operator+(const vec& b); 
	vec& operator-(const vec& b);
	vec& operator*(const vec& b);
	vec& operator/(const vec& b);
	vec& operator=(const vec& b);

	//Elementwise operations vec-scalar
	vec& operator+(const double& b);
	vec& operator-(const double& b);
	vec& operator*(const double& b);
	vec& operator/(const double& b);

	vec transpose();
	double max();
	double min();
	const vec getBlock(int i, int j, int x, int y) const;
	void setBlock(int i, int j, vec v);
	int rows() const;
	int cols() const;
	double norm();
	void resize(int i, int j);

	vec matMult(vec b);
	void setZero();
};
vec sqrt(vec v);
bool sizeCheck(const vec& v1,const vec& v2);

#endif
