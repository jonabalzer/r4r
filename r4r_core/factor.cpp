/*
 * factor.cpp
 *
 *  Created on: Jun 11, 2012
 *      Author: jbalzer
 */

#include "factor.h"
#include <iostream>

extern "C" void dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info);
extern "C" void dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
extern "C" void dpotri_(char* uplo, int* n, double* a, int* lda, int* info);

using namespace std;

namespace R4R {

bool CMatrixFactorization::SVD(const CDenseArray<double>& A, CDenseArray<double>& U, CDenseArray<double>& S, CDenseArray<double>& Vt) {

    // check dimensions
    if(U.NCols()!=A.NRows() || U.NRows()!=A.NRows() || Vt.NCols()!=A.NCols() || Vt.NRows()!=A.NCols()) {

    	cout << "ERROR: Dimension mismatch." << endl;
    	return 1;

    }

	// init
    int m, n, lda, ldu, ldvt, info, lwork;
    m = A.NRows();
    n = A.NCols();
    lda = m;
    ldu = m;
    ldvt = n;
    lwork = -1;

    double wkopt;
    double* work;

    double* a = A.Data().get();
    double* s = S.Data().get();
    double* u = U.Data().get();
    double* vt = Vt.Data().get();

    dgesvd_("A","A",&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,&wkopt,&lwork,&info);

    lwork = (int)wkopt;
    work = (double*)malloc(lwork*sizeof(double));

    dgesvd_("All","All",&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);

    if(info>0) {

    	cout << "ERROR: The algorithm computing SVD failed to converge." << endl;
    	return 1;

    }

	return 0;

}


size_t CMatrixFactorization::Rank(const CDenseArray<double>& A, double tol) {

	CDenseArray<double> U(A.NRows(),A.NRows());
	CDenseArray<double> S(min(A.NRows(),A.NCols()),1);
	CDenseArray<double> Vt(A.NCols(),A.NCols());

	SVD(A,U,S,Vt);

	size_t rank = 0;

	for(size_t i=0; i<S.NRows(); i++) {

		if(S(i,0)>tol)
			rank++;

	}

	return rank;

}

bool CMatrixFactorization::Cholesky(CDenseArray<double>& A) {

   if(A.NCols()!=A.NRows()) {

    	cout << "ERROR: Matrix is not symmetric." << endl;
    	return 1;

    }

   int lda, info, n;
   lda = A.NRows();
   n = lda;
   info = 0;

   double* a = A.Data().get();

   dpotrf_("U",&n,a,&lda,&info);

   if(info>0) {

   	cout << "ERROR: Cholesky decomposition failed." << endl;
   	return 1;

   }

   return 0;

}

bool CMatrixFactorization::InvertSymmetric(CDenseArray<double>& A) {

   if(A.NCols()!=A.NRows()) {

    	cout << "ERROR: Matrix is not symmetric." << endl;
    	return 1;

    }

   Cholesky(A);

   int lda, info, n;
   lda = A.NRows();
   n = lda;
   info = 0;

   double* a = A.Data().get();

   dpotri_("U",&n,a,&lda,&info);

   if(info>0) {

   	cout << "ERROR: Inversion failed." << endl;
   	return 1;

   }

   return 0;

}

}

