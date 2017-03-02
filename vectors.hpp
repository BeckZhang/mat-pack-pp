#include <iostream>
#include <string>
#include <vector>
#include "exceptions.h"

using namespace std;

template<typename TX>
class my_vec {
protected:
    vector<TX> eles;

public:
    my_vec(int n) {
        eles.clear();
        eles.assign(n, 0);
    }

    template<typename X>
    my_vec(my_vec<X> * ob_vec) {
        eles.clear();
        for (int i=0; i<ob_vec->length(); ++i) {
            eles.push_back((TX)(ob_vec->element_access(i+1)));
        }
    }
    
    int length() {
        return eles.size();
    }

    int set_length(int n) {
        if (n > eles.size()) {
            for (int i=eles.size(); i<n; ++i) {
                eles.push_back(0);
            }
            return;
        }
        if (n<eles.size()) {
            eles.erase(eles.begin()+n, eles.end());
        }
    }

    
    virtual void print()=0;
    typename vector<TX>::reference element_access(int n) {
        if (n == -1) {
            n = this->length();
        }

        if (n>this->length() || n<=0) {
            throw InvalidCoordinatesException("index out of range");
        }
        return eles.at(n-1);
    }
    
    typename vector<TX>::reference operator ()(int n) {
        return this->element_access(n);
    }
    
    void append(TX end_ele) {
        this->eles.push_back(end_ele);
    }
    
    void insert(TX new_ele, int pos) {
        if (pos <= 0 || pos > this->length() + 1) {
            throw InvalidCoordinatesException("index out of range");
        }
        this->eles.insert(this->eles.begin() + pos - 1, new_ele);
    }
    
    void remove(int pos) {
        if (pos <= 0 || pos > this->length()) {
            throw InvalidCoordinatesException("index out of range");
        }
        this->eles.erase(this->eles.begin() + pos -1);
    }
};

template<typename TX>
class rowvec;

template<typename TX>
class colvec : public my_vec<TX> {
private:
    colvec<TX> multiply_scalar(TX scal) {
        colvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = scal * this->element_access(i);
        }
        return res;
    }

    colvec<TX> div_scalar(TX scal) {
        colvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = this->element_access(i) / scal;
        }
        return res;
    }
    
    // ============ add, minus, etc ====================
    colvec<TX> add_colvec(colvec<TX> &ob_cvec) {
        if (this->length() != ob_cvec.length()) {
            throw InvalidDimensionsException("sizes of vectors mismatch!");
        }
        colvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res.element_access(i) = this->element_access(i) + ob_cvec.element_access(i);
        }
        return res;
    }
    
    colvec<TX> neg_colvec() {
        colvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = -this->element_access(i);
        }
        return res;
    }
    
public:
    colvec(int n) : my_vec<TX>(n) {;}
    
    template<typename X>
    colvec(const colvec<X> &cvecX) : my_vec<TX>((my_vec<X>*)&cvecX) {;}
    
    void print() {
        cout<<"--"<<endl;
        for (int i=0; i<this->length(); ++i) {
            cout<<this->eles[i]<<endl;
        }
        cout<<"--"<<endl;
    }

    rowvec<TX> T() {
        rowvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = this->element_access(i);
        }
        return res;
    }
    
    //============  OPERATORS ==========================
    colvec<TX> operator + (colvec<TX> &ob_cvec) {
        return this->add_colvec(ob_cvec);
    }
    
    colvec<TX> operator - () {
        return this->neg_colvec();
    }
    
    colvec<TX> operator - (colvec<TX> &ob_cvec) {
        colvec<TX> neg_cvec = ob_cvec.neg_colvec();
        return this->add_colvec(neg_cvec);
    }

    colvec<TX> operator / (TX scal) {
        return this->div_scalar(scal);
    }
    
    //============   FRIENDS  ==========================
    template<typename X>
    friend colvec<X> operator * (X scal, colvec<X> &cvec);
    
};

template<typename TX>
class rowvec : public my_vec<TX> {
private:
    TX right_multiply(colvec<TX> &cvec) {
        if (this->length() != cvec.length()) {
            throw InvalidDimensionsException("sizes of vectors mismatch");
        }
        TX res = 0;
        for (int i = 1; i<=this->length(); ++i) {
            res += this->element_access(i) * cvec(i);
        }
        return res;
    }

    
    rowvec<TX> multiply_scalar(TX scal) {
        rowvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = scal * this->element_access(i);
        }
        return res;
    }

    rowvec<TX> div_scalar(TX scal) {
        rowvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = this->element_access(i) / scal;
        }
        return res;
    }
    
    // ============ add, minus, etc ====================
    rowvec<TX> add_rowvec(rowvec<TX> &ob_rvec) {
        if (this->length() != ob_rvec.length()) {
            throw InvalidDimensionsException("sizes of vectors mismatch!");
        }
        rowvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res.element_access(i) = this->element_access(i) + ob_rvec.element_access(i);
        }
        return res;
    }
    
    rowvec<TX> neg_rowvec() {
        rowvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = -this->element_access(i);
        }
        return res;
    }
public:
    rowvec(int n) : my_vec<TX>(n) {;}
    
    template<typename X>
    rowvec(const rowvec<X> &rvecX) : my_vec<TX>((my_vec<X>*)&rvecX) {;}
    
    void print() {
        cout<<"[";
        for (int i=0; i<this->length(); ++i) {
            if (i==0) {
                cout<<this->eles[0];
            } else {
                cout<<",\t"<<this->eles[i];
            }
        }
        cout<<"]"<<endl;
    }

    colvec<TX> T() {
        colvec<TX> res(this->length());
        for (int i=1; i<=this->length(); ++i) {
            res(i) = this->element_access(i);
        }
        return res;
    }

    //============  OPERATORS  =========================

    TX operator *(colvec<TX> &cvec) {
        return this->right_multiply(cvec);
    }

    rowvec<TX> operator / (TX scal) {
        return this->div_scalar(scal);
    }
    
    rowvec<TX> operator + (rowvec<TX> &ob_rvec) {
        return this->add_colvec(ob_rvec);
    }
    
    rowvec<TX> operator - () {
        return this->neg_rowvec();
    }
    
    rowvec<TX> operator - (rowvec<TX> &ob_rvec) {
        rowvec<TX> neg_rvec = ob_rvec.neg_rowvec();
        return this->add_rowvec(neg_rvec);
    }
    //============   FRIENDS  ==========================
    template<typename X>
    friend rowvec<X> operator * (X scal, rowvec<X> &rvec);
};

template<typename TX>
rowvec<TX> operator * (TX scal, rowvec<TX> &rvec) {
    return rvec.multiply_scalar(scal);
}


template<typename TX>
colvec<TX> operator * (TX scal, colvec<TX> &rvec) {
    return rvec.multiply_scalar(scal);
}

colvec<int> range(int a, int h, int b = -0X3F3F3F3F) {
    if (b == -0X3F3F3F3F) {
        b = h;
        h = 1;
    }
    if (h == 0) {
        throw InvalidCoordinatesException("h can not be zero");
    }
    colvec<int> res(0);
    int i=0;
    while (true) {
        if ( h*(b - (a+h*i)) < 0 ) {
            break;
        }
        res.append(a+h*i);
        i++;
    }
    return res;
}
