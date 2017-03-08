/* This file is part of the MAT-PACK-PP library
 * license  sms, PKU
 * author   Beck Zhang
 * link     https://github.com/BeckZhang/mat-pack-pp
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <iostream>
#include <vector>
#include <string>
#include "vectors.hpp"

#define INF 0X3F3F3F3F
const int END = -1;

using namespace std;

// ==========  INTERFACE  =================
template<typename TX>
class matrix {
private:
    int row_num;
    int col_num;
    vector<vector<TX> > eles;

protected:
    typename vector<TX>::reference element_access(int a, int b);

    matrix<TX> add_matrix(matrix<TX> & mat);
    matrix<TX> minus_matrix(matrix<TX> &mat);
    matrix<TX> neg_matrix();
    matrix<TX> right_multiply_matrix(matrix<TX> &mat);
    matrix<TX> multiply_scalar(TX scal);
    matrix<TX> div_scalar(TX scal);
    matrix<TX> transpose();
    matrix<TX> joint_right(matrix<TX> &mat);
    matrix<TX> joint_right(colvec<TX> &cvec);
    colvec<TX> right_multiply_colvec(colvec<TX> &cvec);

public:
    matrix(int m, int n);
    template<typename X>
        matrix(matrix<X> &mat);
    void print();
    int row_size();
    int col_size();
    void set_row_size(int m);
    void set_col_size(int n);
    void remove_col(int pos);
    void remove_row(int pos);

    matrix<TX> joint_bottom(matrix<TX> mat);
    matrix<TX> joint_bottom(rowvec<TX> rvec);
    matrix<TX> insert_col(colvec<TX> cvec, int pos);
    matrix<TX> insert_row(rowvec<TX> rvec, int pos);

    matrix<TX> submatrix(colvec<int> row_indices, rowvec<int> col_indices);
    matrix<TX> submatrix(int row_start, int row_end, int col_start, int col_end);

    typename vector<TX>::reference operator () (int a, int b);
    matrix<TX> operator +(matrix<TX> &mat);
    matrix<TX> operator -();
    matrix<TX> operator -(matrix<TX> &mat);
    matrix<TX> operator *(matrix<TX> &mat);
    matrix<TX> operator *(TX scal);
    matrix<TX> operator /(TX scal);
    matrix<TX> operator [](matrix<TX> &mat);
    matrix<TX> operator [](colvec<TX> &cvec);
    matrix<TX> T();
	
	colvec<TX> operator * (colvec<TX> &cvec);
	
// ==========  FRIEND FUNCTIONS  ===================
    template<typename X>
        friend matrix<X> operator *(X scal, matrix<X> mat);

    template<typename X>
        friend rowvec<X> operator *(rowvec<X> &rvec, matrix<X> &mat);
};

#endif
