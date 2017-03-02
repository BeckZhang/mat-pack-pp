#include "matrix.cpp"

int main() {
    
    //====  ELEMENT ACCESS & Multiply ====
    matrix<int> mat(3,2);
    colvec<int> cvec(2);
    rowvec<int> rvec(3);

    mat(1,1) = 1;
    mat(1,2) = 2;
    mat(2,1) = 3;
    mat(2,2) = 4;
    mat(3,1) = 5;
    mat(3,2) = 6;


    cvec(1) = 2;
    cvec(2) = 2;

    rvec(1) = 2;
    rvec(2) = 2;
    rvec(3) = 2;

    print("mat = ");
    print( mat );

    print("cvec = ");
    print( cvec );

    print("rvec = ");
    print( rvec );
    print("===========");

    print("mat*cvec = ");
    print( mat*cvec );

    print("rvec*mat = ");
    print( rvec*mat );
    
    //======== INSERT =========
    cvec.insert(100, 3);
    mat = mat[cvec];  // similar with 'mat = [mat, cvec]' in matlab
    print("mat = mat[cvec]");
    print(mat);

    mat = mat.joint_bottom(rvec / 2);
    print("mat = mat.joint_bottom(rvec)");
    print(mat);
    
    //======== HELPERS =======
    print(range(1,10).T());
    print(range(1,-2,-10).T());

    return 0;
}
