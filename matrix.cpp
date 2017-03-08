/* This file is part of the MAT-PACK-PP library
 * license  sms, PKU
 * author   Beck Zhang
 * link     https://github.com/BeckZhang/mat-pack-pp
 */

#include "matrix.h"
#include "exceptions.h"

// ===============  CONTAINS  ======================
template<typename TX>
matrix<TX>::matrix(int m, int n) {
    eles.clear();
    row_num = m;
    col_num = n;
    for (int i=0; i<m; ++i) {
        vector<TX> tmp_vec(n, 0);
        eles.push_back(tmp_vec);
    }
}

template<typename TX>
template<typename X>
matrix<TX>::matrix(matrix<X> &mat) {
    eles.clear();
    row_num = mat.row_size();
    col_num = mat.col_size();
    for (int i=0; i<row_num; ++i) {
        vector<TX> tmp_vec;
        for (int j=0; j<col_num; ++j) {
            tmp_vec.push_back((TX)(mat(i+1, j+1)));
        }
        eles.push_back(tmp_vec);
    }
}

template<typename TX>
void matrix<TX>::print() {
    cout<<"============"<<endl;
    for (int i=0; i<eles.size(); ++i) {
        for (int j=0; j<eles[i].size(); ++j) {
            cout<<eles[i][j]<<",\t";
        }
        cout<<endl;
    }
    cout<<"============"<<endl;
}

template<typename TX>
int matrix<TX>::row_size() {
    row_num = eles.size();
    return row_num;
}

template<typename TX>
int matrix<TX>::col_size() {
    col_num = eles[0].size();
    return col_num;
}

template<typename TX>
void matrix<TX>::set_row_size(int m) {
    if (m>row_num) {
        for (int i=row_num; i<m; ++i) {
            vector<TX> tmp_vec(col_num, 0);
            eles.push_back(tmp_vec);
        }
    } else if (m<row_num) {
        eles.erase(eles.begin()+m, eles.end());
    }
    row_num = m;
}

template<typename TX>
void matrix<TX>::set_col_size(int n) {
    for (int i=0; i<eles.size(); ++i) {
        if (n>col_num) {
            for (int j=col_num; j<n; ++j) {
                eles[i].push_back(0);
            }
        } else if (n<col_num) {
            eles[i].erase(eles[i].begin()+n, eles[i].end());
        }
    }
    col_num = n;
}

template<typename TX>
typename vector<TX>::reference matrix<TX>::element_access(int a, int b) {
    if (a == -1) {
        a = eles.size();
    }
    if (b == -1) {
        b = eles[0].size();
    }

    if (a > row_num || a <= 0) {
        throw InvalidCoordinatesException("Index out of range");
    }
    if (b > col_num || b <= 0) {
        throw InvalidCoordinatesException("Index out of range");
    }
    return eles[a-1].at(b-1);
}

template<typename TX>
typename vector<TX>::reference matrix<TX>::operator () (int a, int b) {
    return this->element_access(a, b);
}

template<typename TX>
matrix<TX> matrix<TX>::add_matrix(matrix<TX> &mat) {
    if (col_num != mat.col_num || row_num != mat.row_num) {
        throw InvalidDimensionsException("sizes of matrix mismatch");
    }
    matrix<TX> res(row_num, col_num);
    for (int i=1; i<=row_num; ++i) {
        for (int j=1; j<=col_num; ++j) {
            res(i,j) = eles[i-1][j-1] + mat(i,j);
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::operator +(matrix<TX> &mat) {
    return this->add_matrix(mat);
}

template<typename TX>
matrix<TX> matrix<TX>::neg_matrix() {
    matrix<TX> res(row_num, col_num);
    for (int i=1; i<=row_num; ++i) {
        for (int j=1; j<=col_num; ++j) {
            res(i,j) = -eles[i-1][j-1];
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::operator -() {
    return this->neg_matrix();
}

template<typename TX>
matrix<TX> matrix<TX>::minus_matrix(matrix<TX> &mat) {
    matrix<TX> neg_mat = mat.neg_matrix();
    return this->add_matrix(neg_mat);
}

template<typename TX>
matrix<TX> matrix<TX>::operator - (matrix<TX> &mat) {
    return this->minus_matrix(mat);
}

template<typename TX>
matrix<TX> matrix<TX>::right_multiply_matrix(matrix<TX> &mat) {
    if (col_num != mat.row_num) {
        throw InvalidDimensionsException("size of matrices dismatch");
    }
    matrix<TX> res(row_num, mat.col_num);
    for (int i=1; i<=row_num; ++i) {
        for (int j=1; j<=mat.col_num; ++j) {
            TX sum_T = 0;
            for (int k=1; k<=col_num; ++k) {
                sum_T = sum_T + (eles[i-1][k-1] * mat(k,j));
            }
            res(i,j) = sum_T;
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::operator * (matrix<TX> &mat) {
    return this->right_multiply_matrix(mat);
}

template<typename TX>
matrix<TX> matrix<TX>::multiply_scalar(TX scal) {
    matrix<TX> res(this->row_size(), this->col_size());
    for (int i=0; i<eles.size(); ++i) {
        for (int j=0; j<eles[0].size(); ++j) {
            res.eles[i][j] = scal * eles[i][j];
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::div_scalar(TX scal) {
    matrix<TX> res(this->row_size(), this->col_size());
    for (int i=0; i<eles.size(); ++i) {
        for (int j=0; j<eles[0].size(); ++j) {
            res.eles[i][j] = eles[i][j] / scal;
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::operator * (TX scal) {
    return this->multiply_scalar(scal);
}

template<typename TX>
matrix<TX> matrix<TX>::operator / (TX scal) {
    return this->div_scalar(scal);
}

template<typename TX>
matrix<TX> matrix<TX>::transpose() {
    matrix<TX> res(this->col_num, this->row_num);
    for (int i=0; i<res.row_num; ++i) {
        for (int j=0; j<res.col_num; ++j) {
            res.eles[i][j] = this->eles[j][i];
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::T() {
    return this->transpose();
}

template<typename TX>
void matrix<TX>::remove_col(int pos) {
    if (pos == -1) {
        pos = this->col_size();
    }
    for (int i=0; i<this->row_size(); ++i) {
        this->eles[i].erase(this->eles[i].begin() + pos -1);
    }
}

template<typename TX>
void matrix<TX>::remove_row(int pos) {
    if (pos == -1) {
        pos = this->row_size();
    }
    this->eles.erase(this->eles.begin() + pos -1);
}

template<typename TX>
matrix<TX> matrix<TX>::joint_right(matrix<TX> &mat) {
    if (this->row_size() != mat.row_size()) {
        throw InvalidDimensionsException("size of matrices mismatch");
    }
    matrix<TX> res(this->row_size(), this->col_size() + mat.col_size());
    for (int i=1; i<=this->row_size(); ++i) {
        for (int j=1; j<=this->col_size(); ++j) {
            res(i,j) = this->element_access(i,j);
        }
        for (int j=1; j<=mat.col_size(); ++j) {
            res(i, j+this->col_size()) = mat(i,j);
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::submatrix(colvec<int> row_indices, rowvec<int> col_indices) {
    //Check Data
    for (int i=1; i<=row_indices.length(); ++i) {
        if (row_indices(i) <= 0 || row_indices(i) > this->row_size()) {
            throw InvalidCoordinatesException("indices out of range");
        }
    }
    for (int j=1; j<=col_indices.length(); ++j) {
        if (col_indices(j) <= 0 || col_indices(j) > this->col_size()) {
            throw InvalidCoordinatesException("indices out of range");
        }
    }
    //return data
    matrix<TX> res(row_indices.length(), col_indices.length());
    for (int i=1; i<=row_indices.length(); ++i) {
        for (int j=1; j<=col_indices.length(); ++j) {
            res(i,j) = this->element_access(row_indices(i), col_indices(j));
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::submatrix(int row_start, int row_end, int col_start, int col_end) {
    if (row_end == -1) {
        row_end = this->row_size();
    }
    if (col_end == -1) {
        col_end = this->col_size();
    }
    if (row_start>row_end || col_start>col_end) {
        throw InvalidCoordinatesException("start index shouldn't be larger than end index"); 
    }
    return this->submatrix(range(row_start, row_end), range(col_start, col_end).T());
}

template<typename TX>
matrix<TX> matrix<TX>::joint_right(colvec<TX> &cvec) {
    if (this->row_size() != cvec.length()) {
        throw InvalidDimensionsException("size of matrix and colvec mismatch");
    }
    matrix<TX> res(this->row_size(), this->col_size() + 1);
    res.eles = this->eles;
    for (int i=1; i<=this->row_size(); ++i) {
        res.eles[i-1].push_back(cvec(i));
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::joint_bottom(matrix<TX> mat) {
    if (this->col_size() != mat.col_size()) {
        throw InvalidDimensionsException("size of matrices mismatch");
    }
    matrix<TX> res(this->row_size() + mat.row_size(), this->col_size());
    for (int i=1; i<=this->row_size(); ++i) {
        for (int j=1; j<=this->col_size(); ++j) {
            res(i,j) = this->element_access(i,j);
        }
    }
    for (int i=1; i<=mat.row_size(); ++i) {
        for (int j=1; j<=this->col_size(); ++j) {
            res(i + this->row_size(), j) = mat(i,j);
        }
    }
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::joint_bottom(rowvec<TX> rvec) {
    if (this->col_size() != rvec.length()) {
        throw InvalidDimensionsException("size of matrix and rowvec mismatch");
    }
    matrix<TX> res(*this);
    vector<TX> tmp_vec;
    for (int i=1; i<=rvec.length(); ++i) {
        tmp_vec.push_back(rvec(i));
    }
    res.eles.push_back(tmp_vec);
    res.row_size();
    res.col_size();
    return res;
}
    
template<typename TX>
matrix<TX> matrix<TX>::operator [] (matrix<TX> &mat) {
    return this->joint_right(mat);
}

template<typename TX>
matrix<TX> matrix<TX>::operator [] (colvec<TX> &cvec) {
    return this->joint_right(cvec);
}

template<typename TX>
matrix<TX> matrix<TX>::insert_col(colvec<TX> cvec, int pos) {
    if (this->row_size() != cvec.length()) {
        throw InvalidDimensionsException("size of matrix and colvec mismatch");
    }

    if (pos == -1) {
        pos = this->col_size() + 1;
    }

    if (pos <= 0 || pos > this->col_size() + 1) {
        throw InvalidCoordinatesException("index out of range");
    }
    matrix<TX> res(*this);
    for (int i=1; i<=this->row_size(); ++i) {
        res.eles[i-1].insert(res.eles[i-1].begin() + pos - 1, cvec(i));
    }
    res.col_size();
    return res;
}

template<typename TX>
matrix<TX> matrix<TX>::insert_row(rowvec<TX> rvec, int pos) {
    if (this->col_size() != rvec.length()) {
        throw InvalidDimensionsException("size of matrix and rowvec mismatch");
    }

    if (pos == -1) {
        pos = this->row_size() + 1;
    }

    if (pos <= 0 || pos > this->row_size() + 1) {
        throw InvalidCoordinatesException("Index out or range");
    }
    matrix res(*this);
    vector<TX> tmp_vec;
    for (int j=1; j<=this->col_size(); ++j) {
        tmp_vec.push_back(rvec(j));
    }
    res.eles.insert(res.eles.begin() + pos -1, tmp_vec);
    res.row_size();
    return res;
}

template<typename TX>
colvec<TX> matrix<TX>::right_multiply_colvec(colvec<TX> &cvec) {
    if (this->col_size() != cvec.length()) {
        throw InvalidDimensionsException("size of matrix and vector mismatch");
    }
    colvec<TX> res(this->row_size());
    for (int i=1; i<=res.length(); ++i) {
        TX sum_TX = 0;
        for (int j=1; j<=this->col_size(); ++j) {
            sum_TX = sum_TX + this->element_access(i,j) * cvec(j);
        }
        res(i) = sum_TX;
    }
    return res;
}

template<typename TX>
colvec<TX> matrix<TX>::operator * (colvec<TX> &cvec) {
	return this->right_multiply_colvec(cvec);
}

// ==========  FRIEND FUNCTIONS  =================
template<typename TX>
matrix<TX> operator * (TX scal, matrix<TX> mat) {
    return mat.multiply_scalar(scal);
}

template<typename TX>
rowvec<TX> operator * (rowvec<TX> &rvec, matrix<TX> &mat) {
    if (rvec.length() != mat.row_size()) {
        throw InvalidDimensionsException("size of colvec and matrix dismatch");
    }
    rowvec<TX> res(mat.col_size());
    for (int j=1; j<=res.length(); ++j) {
        TX sum_TX = 0;
        for (int i=1; i<=mat.row_size(); ++i) {
            sum_TX = sum_TX + rvec(j)*mat(i,j);
        }
        res(j) = sum_TX;
    }
    return res;
}

// ==========  HELPERS  ==========================
inline matrix<int> zeros(int m, int n=-1) {
    if (n == -1) {
        n = m;
    }
    matrix<int> res(m,n);
    return res;
}

inline matrix<int> ones(int m, int n=-1) {
    if (n == -1) {
        n = m;
    }
    matrix<int> res(m,n);
    for (int i=1; i<=m; i++) {
        for (int j=1; j<=n; j++) {
            res(i,j) = 1;
        }
    }
    return res;
}

inline matrix<int> eye(int n) {
    matrix<int> res(n,n);
    for (int i=1; i<=n; ++i) {
        res(i,i) = 1;
    }
    return res;
}

template<typename TX>
inline void print(TX ob) {
    cout<<ob<<endl;
}

template<typename TX>
inline void print(matrix<TX> mat) {
    mat.print();
}

template<typename TX>
inline void print(colvec<TX> cvec) {
    cvec.print();
}

template<typename TX>
inline void print(rowvec<TX> rvec) {
    rvec.print();
}
