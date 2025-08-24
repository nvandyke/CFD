
#ifndef BLOCK_H
#define BLOCK_H

#include <iostream>

class Block {
public:
    Block(int, int);
    Block(int*, int, int);
    Block();
    ~Block();
    Block(const Block&);
    Block& operator=(const Block&);

    inline int& operator()(int x, int y) { return p[x * cols_ + y]; }

    Block& operator+=(const Block&);
    Block& operator-=(const Block&);
    //Block& operator*=(const Block&);
    Block& operator*=(int);
    Block& operator/=(int);

    friend std::ostream& operator<<(std::ostream&, const Block&);
    friend std::istream& operator>>(std::istream&, Block&);

    Block transpose();

    int rows() const;
    int cols() const;
    void print(std::ostream& os = std::cout);
    void setBlock(int i, int j, Block& m);
    const Block getBlock(int i, int j, int x, int y) const;
    int max();
    int size();

private:

    int rows_, cols_;
    int* p;

    void allocSpace();



};

Block operator+(const Block&, const Block&);
Block operator-(const Block&, const Block&);
//Block operator*(const Block&, const Block&);
Block operator*(const Block&, int);
Block operator*(int, const Block&);
Block operator/(const Block&, int);
//Block operator%(Block&, Block&);

#endif
