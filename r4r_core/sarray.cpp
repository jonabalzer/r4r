/*////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2013, Jonathan Balzer
//
// All rights reserved.
//
// This file is part of the R4R library.
//
// The R4R library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The R4R library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with the R4R library. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////*/

#include "sarray.h"
#include "darray.h"

#include <string.h>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include <fstream>


#ifdef HAVE_TBB
#include <tbb/tbb.h>
#endif

using namespace std;

namespace R4R {

template <class T>
CSparseArray<T>::CSparseArray():
	m_nrows(0),
	m_ncols(0),
	m_transpose(false),
    m_data(new spdata) {

}

template <class T>
CSparseArray<T>::CSparseArray(size_t nrows, size_t ncols):
	m_nrows(nrows),
	m_ncols(ncols),
	m_transpose(false),
    m_data(new spdata) {

}

template <class T>
CSparseArray<T>::CSparseArray(const CSparseArray<T>& x):
    m_nrows(x.m_nrows),
    m_ncols(x.m_ncols),
    m_transpose(x.m_transpose),
    m_data(x.m_data) {}

template <class T>
CSparseArray<T>::CSparseArray(size_t nrows, size_t ncols, shared_ptr<spdata> data):
    m_nrows(nrows),
    m_ncols(ncols),
    m_transpose(false),
    m_data(data) {}

template <class T>
CSparseArray<T>& CSparseArray<T>::operator =(const CSparseArray<T>& x) {

    if(this != &x) {

        m_nrows = x.m_nrows;
        m_ncols = x.m_ncols;
        m_transpose = x.m_transpose;
        m_data = x.m_data;
    }

    return *this;

}


template <class T>
CSparseArray<T> CSparseArray<T>::Clone() const {

    // make copy of data
    spdata* data = new spdata;
    *data = *m_data;

    CSparseArray<T> result(m_nrows,m_ncols,shared_ptr<spdata>(data));

    if(m_transpose)
        result.Transpose();

    return result;

}


template <class T>
void CSparseArray<T>::Eye() {

    if(this->Nonzeros()) {

        typename map<size_t,map<size_t,T> >::iterator it;

        for(it=m_data->begin(); it!=m_data->end(); ++it)
            m_data->erase(it);

    }

    for(size_t i=0; i<min(m_nrows,m_ncols); i++)
        this->Set(i,i,1.0);

}

template <class T>
void CSparseArray<T>::Concatenate(const CSparseArray& array, bool direction) {

    // vertical cat
    if(!direction) {

        assert(m_ncols==array.NCols());

        // save old number of rows and resize
        size_t nrowsold = m_nrows;
        this->Resize(m_nrows+array.NRows(),m_ncols);

        typename map<size_t,map<size_t,T> >::const_iterator it;
        typename map<size_t,T>::const_iterator it_col;

        // if both are not transposed, stack row maps
        if(!m_transpose && !array.m_transpose) {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it)
                (*m_data)[nrowsold+it->first] = it->second;

        } // go through rows of input linearly, input col = output row
        else if(m_transpose && !array.m_transpose) {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it) {

                for(it_col=it->second.begin(); it_col!=it->second.end(); ++it_col)
                    (*m_data)[it_col->first][nrowsold+it->first]= it_col->second;

            } // go through input maps (cols) and inject into row maps
        }
        else if(!m_transpose && array.m_transpose) {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it) {

                for(it_col=it->second.begin(); it_col!=it->second.end(); ++it_col)
                    (*m_data)[nrowsold+it_col->first][it->first] = it_col->second;

            }

        } // both transposed,append rows of input to rows of output
        else {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it) {

                for(it_col=it->second.begin(); it_col!=it->second.end(); ++it_col) {
                    (*m_data)[it->first][nrowsold+it_col->first] = it_col->second;

                }


            }

        }

    } else { // horizontal cat

        assert(m_nrows==array.NRows());

        // save old number of colums and resize
        size_t ncolsold = m_ncols;
        this->Resize(m_nrows,m_ncols+array.NCols());

        typename map<size_t,map<size_t,T> >::const_iterator it;
        typename map<size_t,T>::const_iterator it_col;

        // if both are transposed, stack row maps
        if(m_transpose && array.m_transpose) {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it)
                (*m_data)[ncolsold+it->first] = it->second;

        } // go through rows of input linearly, input col = output row
        else if(!m_transpose && array.m_transpose) {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it) {

                for(it_col=it->second.begin(); it_col!=it->second.end(); ++it_col)
                    (*m_data)[it_col->first][ncolsold+it->first]= it_col->second;

            } // go through input maps (cols) and inject into row maps
        }
        else if(m_transpose && !array.m_transpose) {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it) {

                for(it_col=it->second.begin(); it_col!=it->second.end(); ++it_col)
                    (*m_data)[ncolsold+it_col->first][it->first] = it_col->second;

            }

        } // both not transposed, append rows of input to rows of output
        else {

            for(it=array.m_data->begin(); it!=array.m_data->end(); ++it) {

                for(it_col=it->second.begin(); it_col!=it->second.end(); ++it_col) {
                    (*m_data)[it->first][ncolsold+it_col->first] = it_col->second;

                }


            }

        }

    }

}


template<class U>
ostream& operator << (ostream& os, CSparseArray<U>& x) {

	typename map<size_t,map<size_t,U> >::iterator it_row;
	typename map<size_t,U>::iterator it_col;

    it_row = x.m_data->begin();

    while(it_row != x.m_data->end()) {

		it_col = it_row->second.begin();

		while(it_col!=it_row->second.end()) {

			if(x.m_transpose)
				os << "[" << (int)it_col->first << "," << (int)it_row->first << "] " << (float)it_col->second << endl;
			else
				os << "[" << (int)it_row->first << "," << (int)it_col->first << "] " << (float)it_col->second << endl;

			it_col++;

		}

		os << endl;

		it_row++;

	}

	return os;

}

template <class T>
T CSparseArray<T>::Norm2() {

	T sum = 0;

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = m_data->begin();

    while(it_row != m_data->end()) {

		it_col = it_row->second.begin();

		while(it_col!=it_row->second.end()) {

			sum += (it_col->second)*(it_col->second);

			it_col++;

		}

		it_row++;

	}

	return sqrt(sum);

}

template <class T>
void CSparseArray<T>::Transpose() {

    std::swap(m_nrows,m_ncols);
	m_transpose = !m_transpose;

}

template <class T>
T& CSparseArray<T>::operator()(size_t i, size_t j) {

    assert(i<m_nrows && j<m_ncols);

	if(m_transpose)
        std::swap(i,j);

    return (*m_data)[i][j];

}


template <class T>
T CSparseArray<T>::Get(size_t i, size_t j) const {

    // important: do this before the swap
    assert(i<m_nrows && j<m_ncols);

	if(m_transpose)
        std::swap(i,j);

    typename map<size_t,map<size_t,T> >::const_iterator it_row;

    it_row = m_data->find(i);

    if(it_row == m_data->end())
		return (T)0;

    typename map<size_t,T>::const_iterator it_col;

	it_col = it_row->second.find(j);

	if(it_col == it_row->second.end())
		return (T)0;

	return it_col->second;

}

template <class T>
map<size_t,T> CSparseArray<T>::GetRow(size_t i) const {

    assert(i<m_nrows);

	map<size_t,T> row;

	if(!m_transpose) {

        typename map<size_t,map<size_t,T> >::const_iterator loc = m_data->find(i);


        if(loc!=m_data->end())
            row = loc->second;

	}
	else {

        typename map<size_t,map<size_t,T> >::const_iterator it_col = m_data->begin();
        typename map<size_t,T>::const_iterator it_row;

        while(it_col!=m_data->end()) {

			it_row = it_col->second.begin();

			while(it_row!=it_col->second.end()) {


				if(it_col->second.find(i)!=it_col->second.end())
					row.insert(pair<size_t,T>(it_col->first,it_row->second));

				it_row++;

			}

			it_col++;

		}

	}

	return row;

}


template <class T>
void CSparseArray<T>::Set(size_t i, size_t j, T v) {

    // important: do this before swap
    assert(i<m_nrows && j<m_ncols);

	if(m_transpose)
        std::swap(i,j);

	if(v==0){	// if the new value is zero, delete entry

		typename map<size_t,map<size_t,T> >::iterator it_row;
		typename map<size_t,T>::iterator it_col;

        it_row = m_data->find(i);

        if(it_row == m_data->end())
			return;

		it_col = it_row->second.find(j);

		if(it_col == it_row->second.end())
			return;

		it_row->second.erase(it_col);

		if(it_row->second.empty())
            m_data->erase(it_row);

	}
	else
        (*m_data)[i][j] = v;

}

template <class T>
CSparseArray<T> CSparseArray<T>::operator*(const T& scalar) {

	CSparseArray<T> result = *this;

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = result.m_data->begin();

    while(it_row != result.m_data->end()) {

		it_col = it_row->second.begin();

		while(it_col!=it_row->second.end()) {

			it_col->second *= scalar;

			it_col++;

		}

		it_row++;

	}

	return result;

}


template <class T>
CSparseArray<T> CSparseArray<T>::operator+(const CSparseArray& array) {

	assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

	CSparseArray<T> result = array;		// make a copy of the input

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = m_data->begin();

    while(it_row != m_data->end()) {

		it_col = it_row->second.begin();

		while(it_col!=it_row->second.end()) {

			result(it_row->first,it_col->first) += it_col->second;

			it_col++;

		}

		it_row++;

	}

	return result;

}

template <class T>
CSparseArray<T> CSparseArray<T>::operator-(const CSparseArray& array) {

	assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

	CSparseArray<T> result = array;		// make a copy of the input

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = m_data->begin();

    while(it_row != m_data->end()) {

		it_col = it_row->second.begin();

		while(it_col!=it_row->second.end()) {

			result(it_row->first,it_col->first) -= it_col->second;

			it_col++;

		}

		it_row++;

	}

	return result;

}

template <class T>
T CSparseArray<T>::InnerProduct(const CSparseArray<T>& x, const CSparseArray<T>& y) {

	assert(x.m_nrows==y.m_nrows && x.m_ncols==y.m_ncols);

	T sum = 0;

    typename map<size_t,map<size_t,T> >::const_iterator it_row;
    typename map<size_t,T>::const_iterator it_col;

    it_row = x.m_data->begin();

    while(it_row!=x.m_data->end()) {

		map<size_t,T> rowx = x.GetRow(it_row->first);
		map<size_t,T> rowy = y.GetRow(it_row->first);

		for(it_col=rowx.begin(); it_col!=rowx.end(); it_col++)
			sum+= it_col->second*rowy[it_col->first];

		it_row++;

	}

	return sum;

}

template <class T>
template<class Matrix> Matrix CSparseArray<T>::operator*(const Matrix& array) const {

//#ifdef _OPENMP
//    assert(m_ncols==array.NRows());

//    Matrix result(m_nrows,array.NCols());;

//    if(!m_transpose) {

//        typename map<size_t,T>::iterator it_col;

//#pragma omp parallel for private(it_col), shared(result)
//        for(size_t i=0; i<m_nrows; i++) {

//            for(size_t j=0; j<array.NCols(); j++) {

//                if(m_data[i].size()>0) {

//                    T sum = 0;

//                    for(it_col=m_data[i].begin(); it_col!=m_data[i].end(); it_col++)
//                        sum+= (it_col->second)*array.Get(it_col->first,j);

//                    if(sum!=0)
//                        result(i,j) = sum;

//                }

//            }

//        }

//    }
//    else {

//        typename map<size_t,T>::iterator it_row;

//        for(size_t j=0; j<array.NCols(); j++) {

//#pragma omp parallel for private(it_row), shared(result)
//            for(size_t i=0; i<m_ncols; i++) {


//                if(m_data[i].size()>0) {

//                    for(it_row=m_data[i].begin(); it_row!=m_data[i].end(); it_row++) {




//#pragma omp critical
//                        {
//                        T val = it_row->second*array.Get(i,j);

//                        // check this to avoid fill-in
//                        if(val!=0)
//                            result(it_row->first,j) += val;

//                        }

//                    }

//                }

//            }

//        }

//    }

//    return result;
// #else


	assert(m_ncols==array.NRows());

	Matrix result(m_nrows,array.NCols());;

	if(!m_transpose) {

        typename map<size_t,map<size_t,T> >::const_iterator it_row;
        typename map<size_t,T>::const_iterator it_col;

        for(it_row=m_data->begin(); it_row!=m_data->end(); it_row++) {

			for(size_t j=0; j<array.NCols(); j++) {

				T sum = 0;

				for(it_col=it_row->second.begin(); it_col!=it_row->second.end(); it_col++)
					sum+= (it_col->second)*array.Get(it_col->first,j);

				if(sum!=0)
					result(it_row->first,j) = sum;

			}

		}

	}
	else {

        typename map<size_t,map<size_t,T> >::const_iterator it_col;
        typename map<size_t,T>::const_iterator it_row;

		for(size_t j=0; j<array.NCols(); j++) {

                for(it_col=m_data->begin(); it_col!=m_data->end(); it_col++) {

				for(it_row=it_col->second.begin(); it_row!=it_col->second.end(); it_row++) {

					T val = it_row->second*array.Get(it_col->first,j);

					if(val!=0)
						result(it_row->first,j) += val;

				}

			}

		}

	}

	return result;
//#endif

}

//template <class T>
// CDenseVector<T> CSparseArray<T>::operator*(const CDenseVector<T>& vector) {

//#ifdef _OPENMP
//    assert(m_ncols==vector.NRows());

//    CDenseVector<T> result = CDenseVector<T>(m_nrows);

//    if(!m_transpose) {

//        typename map<size_t,T>::iterator it_col;

//#pragma omp parallel for private(it_col), shared(result,vector)
//        for(size_t i=0; i<m_nrows; i++) {

//            if(m_data[i].size()>0) {

//                T sum = 0;

//                for(it_col=m_data[i].begin(); it_col!=m_data[i].end(); it_col++)
//                    sum+= (it_col->second)*vector.Get(it_col->first);

//                result(i) = sum;

//            }

//        }

//    }
//    else {

//        typename map<size_t,T>::iterator it_row;

//#pragma omp parallel for private(it_row), shared(result,vector)
//        for(size_t i=0; i<m_ncols; i++) {

//            if(m_data[i].size()>0) {

//                for(it_row=m_data[i].begin(); it_row!=m_data[i].end(); it_row++) {

//#pragma omp critical
//                    result(it_row->first) += it_row->second*vector.Get(i);				// protect this

//                }

//            }

//        }

//    }

//    return result;
//#else

//    assert(m_ncols==vector.NRows());

//	CDenseVector<T> result = CDenseVector<T>(m_nrows);

//	if(!m_transpose) {

//		typename map<size_t,map<size_t,T> >::iterator it_row;
//		typename map<size_t,T>::iterator it_col;

//		for(it_row=m_data.begin(); it_row!=m_data.end(); it_row++) {

//			T sum = 0;

//			for(it_col=it_row->second.begin(); it_col!=it_row->second.end(); it_col++)
//				sum+= (it_col->second)*vector.Get(it_col->first);

//			result(it_row->first) = sum;

//		}

//	}
//	else {

//		typename map<size_t,map<size_t,T> >::iterator it_col;
//		typename map<size_t,T>::iterator it_row;

//		for(it_col=m_data.begin(); it_col!=m_data.end(); it_col++) {

//			for(it_row=it_col->second.begin(); it_row!=it_col->second.end(); it_row++)
//				result(it_row->first) += it_row->second*vector.Get(it_col->first);

//		}

//	}

//	return result;
//#endif

//}

template <class T>
void CSparseArray<T>::Scale(T scalar) {

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = m_data->begin();

    while(it_row != m_data->end()) {

		it_col = it_row->second.begin();

		while(it_col!=it_row->second.end()) {

			it_col->second *= scalar;

			it_col++;

		}

		it_row++;

	}

}

template <class T>
void CSparseArray<T>::ScaleRow(size_t i, T scalar) {

    typename map<size_t,map<size_t,T> >::iterator it_row;

    it_row = m_data->find(i);

    if(it_row!=m_data->end()) {

        typename map<size_t,T>::iterator it_col;

        for(it_col=it_row->second.begin(); it_col !=it_row->second.end(); ++it_col)
            it_col->second *= scalar;

    }

}

template <class T>
CSparseArray<T> CSparseArray<T>::Transpose(const CSparseArray<T>& array) {

    // shallow copy
	CSparseArray<T> result = array;

	result.Transpose();

	return result;

}

template <class T>
CSparseArray<T> CSparseArray<T>::Square(const CSparseArray<T> &array) {


	CSparseArray<T> result(array.m_nrows,array.m_nrows);

    typename map<size_t,map<size_t,T> >::const_iterator it_row1;
    typename map<size_t,map<size_t,T> >::const_iterator it_row2;
    typename map<size_t,T>::const_iterator it_el;

	// multiply each row with itself if there is enough overlap
    for(it_row1 = array.m_data->begin(); it_row1 != array.m_data->end(); it_row1++) {

        for(it_row2 = array.m_data->begin(); it_row2 != array.m_data->end(); it_row2++) {

			// check if there is overlap
			if(it_row1->second.rbegin()->first >= it_row2->second.begin()->first) {

				T sum = 0;

				for(it_el = it_row1->second.begin(); it_el != it_row1->second.end(); it_el++)
					sum += it_el->second*array.Get(it_row2->first,it_el->first);   // no: maybe access over it_row2 but make sure to avoid fill-in

				// check this to avoid fill-in
				if(sum!=0)
					result(it_row1->first,it_row2->first) = sum;

			}

		}

	}

	return result;

}




template <class T>
size_t CSparseArray<T>::Nonzeros() {

	typename map<size_t,map<size_t,T> >::iterator it_row;

	size_t nnz = 0;

    it_row = m_data->begin();

    while(it_row!=m_data->end()){

		nnz += it_row->second.size();

		it_row++;

	}

	return nnz;

}

template <class T>
bool CSparseArray<T>::Symmetric() {

	bool result = true;

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = m_data->begin();

    while(it_row!=m_data->end()){

		it_col = it_row->second.begin();

		while(it_col->first<it_row->first) {

			if(it_col->second!=Get(it_col->first,it_row->first)) {

				return false;

			}

			it_col++;

		}

		it_row++;

	}

	return result;

}


template <class T>
void CSparseArray<T>::DeleteRow(size_t i) {

	typename map<size_t,map<size_t,T> >::iterator it_row;

    it_row = m_data->find(i);

    if(it_row!=m_data->end())
        m_data->erase(it_row);

}

template <class T>
void CSparseArray<T>::DeleteColumn(size_t j) {

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = m_data->begin();

    while(it_row!=m_data->end()){

		it_col = it_row->second.find(j);

		if(it_col!=it_row->second.end())
			it_row->second.erase(it_col);

		it_row++;

	}

}

template <class T>
void CSparseArray<T>::GetCSR(vector<size_t>& nz, vector<size_t>& j, vector<T>& v, bool ibase) {

    typename map<size_t,map<size_t,T> >::iterator it_row = m_data->begin();
	typename map<size_t,T>::iterator it_col;

	size_t nnz = ibase;

    while(it_row!=m_data->end()){

		nz.push_back(nnz);

		for (it_col=it_row->second.begin(); it_col!=it_row->second.end(); it_col++) {

			j.push_back(it_col->first+ibase);

			v.push_back((T)(it_col->second));

			nnz++;

		}

		it_row++;

	}

	nz.push_back(nnz);

}

template <class T>
void CSparseArray<T>::GetCOO(std::vector<size_t>& i, std::vector<size_t>& j, std::vector<T>& v, bool ibase) {

    typename map<size_t,map<size_t,T> >::iterator it_row = m_data->begin();
	typename map<size_t,T>::iterator it_col;

    while(it_row!=m_data->end()){

		for (it_col=it_row->second.begin(); it_col!=it_row->second.end(); it_col++) {

			i.push_back(it_row->first+ibase);
			j.push_back(it_col->first+ibase);
			v.push_back((T)(it_col->second));

		}

		it_row++;

	}

}

template <class T>
bool CSparseArray<T>::WriteToFile(const char* filename) {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file.\n";
		return 1;

	 }

	// write header
	out << "\%\%MatrixMarket matrix coordinate real general" << endl;
	out << m_nrows << " " << m_ncols << " " << Nonzeros() << endl;

	typename map<size_t,map<size_t,T> >::iterator it_row;
	typename map<size_t,T>::iterator it_col;

    it_row = m_data->begin();

    while(it_row != m_data->end()) {

		it_col = it_row->second.begin();

		while(it_col!=it_row->second.end()) {

			if(m_transpose)
				out << (int)it_col->first << " " << (int)it_row->first << " " << (float)it_col->second << endl;
			else
				out << (int)it_row->first << " " << (int)it_col->first << " " << (float)it_col->second << endl;

			it_col++;

		}

		it_row++;

	}

	out.close();

	return 0;

}

template <class T>
void CSparseArray<T>::Rand(size_t nnz) {

	for(size_t k=0; k<nnz; k++) {

		size_t i = (size_t)floor(((double)rand()/(double)RAND_MAX)*(double)m_nrows);
		size_t j = (size_t)floor(((double)rand()/(double)RAND_MAX)*(double)m_ncols);
		this->operator ()(i,j) = (T)((T)rand()/(T)RAND_MAX);

	}

}

template class CSparseArray<float>;
template class CSparseArray<double>;
template class CSparseArray<int>;
template CDenseArray<double> CSparseArray<double>::operator *(const CDenseArray<double>& x) const;
template CDenseVector<double> CSparseArray<double>::operator *(const CDenseVector<double>& x) const;
template CSparseArray<double> CSparseArray<double>::operator *(const CSparseArray<double>& x) const;
template CDenseArray<float> CSparseArray<float>::operator *(const CDenseArray<float>& x) const;
template CDenseVector<float> CSparseArray<float>::operator *(const CDenseVector<float>& x) const;
template CSparseArray<float> CSparseArray<float>::operator *(const CSparseArray<float>& x) const;
template ostream& operator<< (ostream& os, CSparseArray<double>& x);
template ostream& operator<< (ostream& os, CSparseArray<float>& x);

template <class T>
CSparseBandedArray<T>::CSparseBandedArray():
	m_nrows(0),
	m_ncols(0),
	m_transpose(false),
	m_data() {

}

template <class T>
CSparseBandedArray<T>::CSparseBandedArray(size_t nrows, size_t ncols):
	m_nrows(nrows),
	m_ncols(ncols),
	m_transpose(false),
	m_data() {

}

template <class T>
CSparseBandedArray<T> CSparseBandedArray<T>::Transpose(const CSparseBandedArray<T>& array) {

	CSparseBandedArray<T> result = array;

	result.Transpose();

	return result;

}

template <class T>
void CSparseBandedArray<T>::Delete(size_t i, size_t j) {

    assert(i<m_nrows && j<m_ncols);

	if(m_transpose)
        std::swap(i,j);

	int b = Band(i,j);
	size_t d = BandIndex(i,j);

	// delete the value if it exists
	m_data[b].erase(d);

	// remove the band, if it is now empty
	if(m_data[b].empty())
		m_data.erase(b);

}

template <class T>
void CSparseBandedArray<T>::DeleteRow(size_t i) {

    assert(i<m_nrows);

	for(size_t j=0; j<m_ncols; j++)
		Delete(i,j);

}

template <class T>
void CSparseBandedArray<T>::DeleteCol(size_t j) {

    assert(j<m_ncols);

	for(size_t i=0; i<m_nrows; i++)
		Delete(i,j);

}

template<class U>
ostream& operator << (ostream& os, const CSparseBandedArray<U>& x) {

    typename map<int,map<size_t,U> >::const_iterator it_b;
    typename map<size_t,U>::const_iterator it_d;

    it_b = x.m_data.begin();

    while(it_b != x.m_data.end()) {

        it_d = it_b->second.begin();

        while(it_d!=it_b->second.end()) {

            if(x.m_transpose)
                os << "[ " << x.Row(-it_b->first,it_d->first) << " " << x.Col(-it_b->first,it_d->first) << " " << it_d->second << " ]" << endl;
            else
                os << "[ " << x.Row(it_b->first,it_d->first) << " " <<  x.Col(it_b->first,it_d->first) << " " << it_d->second << " ]" << endl;

            it_d++;

        }

        os << endl;

        it_b++;

    }

    return os;

}

template <class T>
T CSparseBandedArray<T>::Get(size_t i, size_t j) {

	if(m_transpose)
        std::swap(i,j);

    assert(i<m_nrows && j<m_ncols);

	// convert index
	int b = Band(i,j);
	size_t d = BandIndex(i,j);

	typename map<int,map<size_t,T> >::iterator it_b;

	it_b = m_data.find(b);

	if(it_b == m_data.end())
		return (T)0;

	typename map<size_t,T>::iterator it_d;

	it_d = it_b->second.find(d);

	if(it_d == it_b->second.end())
		return (T)0;

	return it_d->second;

}

template <class T>
map<int,map<size_t,T> > CSparseBandedArray<T>::GetBands(int lower, int upper) {

	assert(lower<=upper);

	if(m_transpose) {

		int temp = lower;
		lower = -upper;
		upper = -temp;

	}

	map<int,map<size_t,T> > result;

	typename map<int,map<size_t,T> >::iterator it_b = m_data.find(lower);

	if(it_b==m_data.end())
		it_b = m_data.begin();

	while(it_b->first <=upper && it_b!=m_data.end()) {

		if(it_b->first>=lower) {

			if(m_transpose)
				result.insert(pair<int,std::map<size_t,T> >(-it_b->first,it_b->second));
			else
				result.insert(pair<int,std::map<size_t,T> >(it_b->first,it_b->second));

		}

		it_b++;
	}

	return result;

}


template <class T>
void CSparseBandedArray<T>::Set(size_t i, size_t j, T v) {

    assert(i<m_nrows && j<m_ncols);

	if(m_transpose)
        std::swap(i,j);

	// convert index
	int b = Band(i,j);
	size_t d = BandIndex(i,j);

    if(v==0) {	// if the new value is zero, delete entry

		typename map<int,map<size_t,T> >::iterator it_b;
		typename map<size_t,T>::iterator it_d;

		it_b = m_data.find(b);

		if(it_b == m_data.end())
			return;

		it_d = it_b->second.find(d);

		if(it_d == it_b->second.end())
			return;

		it_b->second.erase(it_d);

		if(it_b->second.empty())
			m_data.erase(it_b);

	}
	else
		m_data[b][d] = v;

}



template <class T>
T& CSparseBandedArray<T>::operator()(size_t i, size_t j) {

    assert(i<m_nrows && j<m_ncols);

	if(m_transpose)
        std::swap(i,j);

	int b = Band(i,j);
	size_t d = BandIndex(i,j);

	return m_data[b][d];

}

template <class T>
T CSparseBandedArray<T>::Norm2() {

	T sum = 0;

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	it_b = m_data.begin();

	while(it_b != m_data.end()) {

		it_d = it_b->second.begin();

		while(it_d!=it_b->second.end()) {

			sum += (it_d->second)*(it_d->second);

			it_d++;

		}

		it_b++;

	}

	return sqrt(sum);

}

template <class T>
void CSparseBandedArray<T>::Transpose() {

    std::swap(m_nrows,m_ncols);
	m_transpose = !m_transpose;

}

template <class T>
CSparseBandedArray<T> CSparseBandedArray<T>::operator+(const CSparseBandedArray& array) {

	assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

	CSparseBandedArray<T> result = array;		// make a copy of the input

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	it_b = m_data.begin();

	while(it_b != m_data.end()) {

		it_d = it_b->second.begin();

		while(it_d!=it_b->second.end()) {

			if(result.m_transpose)
				result.m_data[-it_b->first][it_d->first] += it_d->second;
			else
				result.m_data[it_b->first][it_d->first] += it_d->second;

			it_d++;

		}

		it_b++;

	}

	return result;

}


template <class T>
CSparseBandedArray<T> CSparseBandedArray<T>::operator-(const CSparseBandedArray& array) {

	assert(m_nrows==array.m_nrows && m_ncols==array.m_ncols);

	CSparseBandedArray<T> result = array;		// make a copy of the input

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	it_b = m_data.begin();

	while(it_b != m_data.end()) {

		it_d = it_b->second.begin();

		while(it_d!=it_b->second.end()) {

			if(result.m_transpose)
				result.m_data[-it_b->first][it_d->first] -= it_d->second;
			else
				result.m_data[it_b->first][it_d->first] -= it_d->second;

			it_d++;

		}

		it_b++;

	}

	return result;

}

template <class T>
CDenseArray<T> CSparseBandedArray<T>::operator*(const CDenseArray<T>& array) {

	assert(m_ncols==array.NRows());

	CDenseArray<T> result = CDenseArray<T>(m_nrows,array.NCols());;

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	// for all columns
	for(size_t k=0; k<array.NCols(); k++) {

		for(it_b=m_data.begin(); it_b!=m_data.end(); it_b++) {

			for(it_d=it_b->second.begin(); it_d!=it_b->second.end(); it_d++) {

				size_t i = Row(it_b->first,it_d->first);
				size_t j = Col(it_b->first,it_d->first);

				if(m_transpose)
                    std::swap(i,j);

				result(i,k) += it_d->second*array.Get(j,k);

			}

		}

	}

	return result;

}

template <class T>
CSparseBandedArray<T> CSparseBandedArray<T>::operator*(CSparseBandedArray<T>& array) {

    assert(m_ncols==array.NRows());

	CSparseBandedArray<T> result = CSparseBandedArray<T>(m_nrows,array.NCols());;

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	// for all columns
	for(size_t k=0; k<array.NCols(); k++) {

		for(it_b=m_data.begin(); it_b!=m_data.end(); it_b++) {

			for(it_d=it_b->second.begin(); it_d!=it_b->second.end(); it_d++) {

				size_t i = Row(it_b->first,it_d->first);
				size_t j = Col(it_b->first,it_d->first);

				if(m_transpose)
                    std::swap(i,j);

				T temp = it_d->second*array.Get(j,k);

				if(temp!=0)
					result(i,k) += temp;

			}

		}

	}

	return result;

}


template <class T>
void CSparseBandedArray<T>::Scale(T scalar) {

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	it_b = m_data.begin();

	while(it_b != m_data.end()) {

		it_d = it_b->second.begin();

		while(it_d!=it_b->second.end()) {

			it_d->second *= scalar;

			it_d++;

		}

		it_b++;

	}

}

template <class T>
void CSparseBandedArray<T>::AddDiagonal(const T scalar) {

	for(size_t i=0; i<min(m_nrows,m_ncols); i++)
		this->operator ()(i,i) += scalar;

}

template <class T>
void CSparseBandedArray<T>::ScaleDiagonal(T scalar) {

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	it_b = m_data.find(0);

	if(it_b==m_data.end())
		return;

	it_d = it_b->second.begin();

	while(it_d!=it_b->second.end()) {

		it_d->second *= scalar;

		it_d++;

	}

}


template <class T>
size_t CSparseBandedArray<T>::Nonzeros() {

	typename map<int,map<size_t,T> >::iterator it_b;

	size_t nnz = 0;

	it_b = m_data.begin();

	while(it_b != m_data.end()) {

		nnz += it_b->second.size();

		it_b++;

	}

	return nnz;

}

template <class T>
bool CSparseBandedArray<T>::SaveToFile(const char* filename) {

	ofstream out(filename);

	if(!out) {

		cout << "ERROR: Could not open file.\n";
		return 1;

	 }

	// write header
	out << "\%\%MatrixMarket matrix coordinate real general" << endl;
	out << m_nrows << " " << m_ncols << " " << Nonzeros() << endl;

	// write data
	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	it_b = m_data.begin();

	while(it_b != m_data.end()) {

		it_d = it_b->second.begin();

		while(it_d!=it_b->second.end()) {

			size_t i = Row(it_b->first,it_d->first);
			size_t j = Col(it_b->first,it_d->first);

			if(m_transpose)
                std::swap(i,j);

			out << i << " " << j << " " << it_d->second;

			it_d++;

		}

		it_b++;

	}

	out.close();

	return 0;

}

template class CSparseBandedArray<float>;
template class CSparseBandedArray<double>;
template ostream& operator<< (ostream& os, const CSparseBandedArray<double>& x);
template ostream& operator<< (ostream& os, const CSparseBandedArray<float>& x);

template <class T>
CSparseDiagonalArray<T>::CSparseDiagonalArray():
	CSparseBandedArray<T>() {


}

template <class T>
CSparseDiagonalArray<T>::CSparseDiagonalArray(size_t nrows, size_t ncols):
	CSparseBandedArray<T>(min(nrows,ncols),min(nrows,ncols)) {

}

template <class T>
T CSparseDiagonalArray<T>::Get(size_t i, size_t j) {

    return j==i ? CSparseBandedArray<T>::Get(i,j) : 0;

}

template <class T>
void CSparseDiagonalArray<T>::Set(size_t i, size_t j, T v) {

    if(j==i)
        CSparseBandedArray<T>::Set(i,j,v);

}

template <class T>
T& CSparseDiagonalArray<T>::operator()(size_t i, size_t j) {

	return CSparseBandedArray<T>::operator()(i,j);

}

template <class T>
void CSparseDiagonalArray<T>::Solve(CDenseArray<T>& x, const CDenseArray<T>& b) const {

    assert(b.NRows()==m_nrows && x.NRows()==m_nrows && b.NCols()==x.NCols());

    typename map<int,map<size_t,T> >::const_iterator it_b;
    typename map<size_t,T>::const_iterator it_d;

	it_b = m_data.find(0);

	if(it_b->second.size()!= b.NRows()) {

		cout << "ERROR: Matrix is rank-deficient." << endl;
		return;

	}

	it_d = it_b->second.begin();

	while(it_d!=it_b->second.end()) {

        for(size_t j=0; j<b.NCols(); j++)
            x(it_d->first,j)=b.Get(it_d->first,j)/it_d->second;

        it_d++;

	}

}


template <class T>
void CSparseDiagonalArray<T>::Invert() {

	typename map<int,map<size_t,T> >::iterator it_b;
	typename map<size_t,T>::iterator it_d;

	it_b = m_data.find(0);

	if(it_b->second.size()!=m_nrows) {

        cerr << "ERROR: Matrix is rank-deficient." << endl;
		return;

	}

	it_d = it_b->second.begin();

	while(it_d!=it_b->second.end()) {

		it_d->second = 1/it_d->second;

		it_d++;

	}

}


template class CSparseDiagonalArray<float>;
template class CSparseDiagonalArray<double>;

template <class T>
CSparseUpperTriangularArray<T>::CSparseUpperTriangularArray():
	CSparseBandedArray<T>() {

}

template <class T>
CSparseUpperTriangularArray<T>::CSparseUpperTriangularArray(size_t nrows, size_t ncols):
	CSparseBandedArray<T>(min(nrows,ncols),min(nrows,ncols)) {

}

template <class T>
T CSparseUpperTriangularArray<T>::Get(size_t i, size_t j) {

	return j>=i ? CSparseBandedArray<T>::Get(i,j) : 0;

}


template <class T>
void CSparseUpperTriangularArray<T>::Set(size_t i, size_t j, T v) {

    if(j>=i)
        CSparseBandedArray<T>::Set(i,j,v);

}

template <class T>
T& CSparseUpperTriangularArray<T>::operator()(size_t i, size_t j) {

	return CSparseBandedArray<T>::operator()(i,j);

}

template <class T>
void CSparseUpperTriangularArray<T>::Solve(CDenseArray<T>& x, const CDenseArray<T>& b) const {

    assert(m_nrows>=m_ncols && b.NRows()==m_nrows && x.NRows()==m_ncols && b.NCols()==x.NCols());

    typename map<int,map<size_t,T> >::const_iterator it_diagonal = m_data.find(0);

	if(it_diagonal->second.size()!= b.NRows()) {

        cerr << "ERROR: Matrix is rank-deficient." << endl;
		return;

	}

    map<int,typename map<size_t,T>::const_reverse_iterator> bands;
    typename map<int,map<size_t,T> >::const_iterator it_b;

	// increment to start adding everyting above diagonal
	it_diagonal++;

	// init iterators to the beginning of all bands
	for(it_b=it_diagonal; it_b!=m_data.end();it_b++)
        bands.insert(pair<int,typename map<size_t,T>::const_reverse_iterator >(it_b->first,it_b->second.rbegin()));

	// set iterator back to diagonal
	it_diagonal--;

    typename map<size_t,T>::const_reverse_iterator it_d = it_diagonal->second.rbegin();
    typename map<int,typename map<size_t,T>::const_reverse_iterator >::iterator it_bands;

	// one accumulator for each column of b
	T* sums = new T[b.NCols()];

	// iterate through main diagonal
	while(it_d!=it_diagonal->second.rend()) {

		for(size_t bc=0; bc<b.NCols(); bc++)
			sums[bc] = b.Get(it_d->first,bc);

		for(it_bands=bands.begin(); it_bands!=bands.end(); it_bands++) {

			// go through all band iterators and compute their row/column representation
			size_t i = CSparseBandedArray<T>::Row(it_bands->first,it_bands->second->first);
			size_t j = CSparseBandedArray<T>::Col(it_bands->first,it_bands->second->first);


			// check whether the band contributes to the current row
			if(i==it_d->first && j>it_d->first) {

				for(size_t bc=0; bc<b.NCols(); bc++)
					sums[bc] -= it_bands->second->second*x(j,bc);

				// increment the iterator in the band
				it_bands->second++;

			}

		}

		for(size_t bc=0; bc<b.NCols(); bc++)
			x(it_d->first,bc) = sums[bc]/it_d->second;


		it_d++;

	}

	delete [] sums;

}


template class CSparseUpperTriangularArray<float>;
template class CSparseUpperTriangularArray<double>;


template <class T>
CSparseLowerTriangularArray<T>::CSparseLowerTriangularArray():
	CSparseBandedArray<T>() {

}

template <class T>
CSparseLowerTriangularArray<T>::CSparseLowerTriangularArray(size_t nrows, size_t ncols):
	CSparseBandedArray<T>(min(nrows,ncols),min(nrows,ncols)) {

}

template <class T>
T CSparseLowerTriangularArray<T>::Get(size_t i, size_t j) {

	return j<=i ? CSparseBandedArray<T>::Get(i,j) : 0;

}


template <class T>
void CSparseLowerTriangularArray<T>::Set(size_t i, size_t j, T v) {

    if(j<=i)
        CSparseBandedArray<T>::Set(i,j,v);

}

template <class T>
T& CSparseLowerTriangularArray<T>::operator()(size_t i, size_t j) {

     // just work with diagonal values, don't worry about memory
    //assert(j<=i);

	return CSparseBandedArray<T>::operator()(i,j);

}

template <class T>
void CSparseLowerTriangularArray<T>::Solve(CDenseArray<T>& x, const CDenseArray<T>& b) const {

    assert(m_nrows>=m_ncols && b.NRows()==m_nrows && x.NRows()==m_ncols && b.NCols()==x.NCols());

    typename map<int,map<size_t,T> >::const_iterator it_diagonal = m_data.find(0);

	if(it_diagonal->second.size()!= b.NRows()) {

        cerr << "ERROR: Matrix is rank-deficient." << endl;
		return;

	}

    map<int,typename map<size_t,T>::const_iterator> bands;
    typename map<int,map<size_t,T> >::const_iterator it_b;

	// init iterators to the beginning of all bands
	for(it_b=m_data.begin(); it_b!=m_data.end();it_b++)
        bands.insert(pair<int,typename map<size_t,T>::const_iterator >(it_b->first,it_b->second.begin()));

    typename map<size_t,T>::const_iterator it_d;
    typename map<int,typename map<size_t,T>::const_iterator >::iterator it_bands;

	it_d = it_diagonal->second.begin();

	// one accumulator for each column of b
	T* sums = new T[b.NCols()];

	// iterate through main diagonal
	while(it_d!=it_diagonal->second.end()) {

		for(size_t bc=0; bc<b.NCols(); bc++)
			sums[bc] = b.Get(it_d->first,bc);

		for(it_bands=bands.begin(); it_bands!=bands.end(); it_bands++) {

			// go through all band iterators and compute their row/column representation
			size_t i = CSparseBandedArray<T>::Row(it_bands->first,it_bands->second->first);
			size_t j = CSparseBandedArray<T>::Col(it_bands->first,it_bands->second->first);

			// check whether the band contributes to the current row
			if(i==it_d->first && j<it_d->first) {

				for(size_t bc=0; bc<b.NCols(); bc++)
					sums[bc] -= it_bands->second->second*x(j,bc);

				// increment the iterator in the band
				it_bands->second++;

			}

		}

		for(size_t bc=0; bc<b.NCols(); bc++) {
			x(it_d->first,bc) = sums[bc]/it_d->second;
			//cout << sums[bc]/it_d->second << endl;
		}


		it_d++;

	}

	delete [] sums;

}

template class CSparseLowerTriangularArray<float>;
template class CSparseLowerTriangularArray<double>;

}


