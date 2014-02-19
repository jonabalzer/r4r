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

#ifndef R4RRBUFFER_H
#define R4RRBUFFER_H

#include <stdlib.h>
#include <vector>
#include <iostream>

namespace R4R {

/*! \brief a simple array-based ring buffer implementation
 *
 * The interface is made STL-compliant (i.e., similar to vectors and lists)
 * to make the container usable as template parameter.
 *
 *
 */
template<class T,class Allocator = std::allocator<T> >
class CRingBuffer {

public:

    class iterator {

    public:

      iterator():m_i(0),m_n(0),m_data(nullptr) {}

      iterator(int i, size_t n, T* data):m_i(i),m_n(n),m_data(data) {}

      void operator++() { m_i = ++m_i%m_n; }

      iterator operator+(int i) const { return iterator((m_i+i)%m_n,m_n,m_data); }

      void operator+=(int i) { m_i=(m_i+i)%m_n; }

      T& operator*() { return m_data[m_i]; }

      T* operator->() { return &(operator*()); }

      bool operator!=(const iterator& it) const { return m_i!=it.m_i; }

    private:

        int m_i;
        size_t m_n;
        T* m_data;

    };

    class const_iterator {

    public:

      const_iterator():m_i(0),m_n(0),m_data(nullptr) {}

      const_iterator(int i, size_t n, const T* data):m_i(i),m_n(n),m_data(data) {}

      void operator++() { m_i = ++m_i%m_n; }

      const_iterator operator+(int i) const { return const_iterator((m_i+i)%m_n,m_n,m_data); }

      void operator+=(int i) { m_i=(m_i+i)%m_n; }

      const T* operator->() const { return &(operator*()); }

      const T& operator*() const { return m_data[m_i]; }

      bool operator!=(const const_iterator& it) const { return m_i!=it.m_i; }

    private:

        int m_i;
        size_t m_n;
        const T* m_data;

    };

    class reverse_iterator {

    public:

      reverse_iterator():m_i(0),m_n(0),m_data(nullptr) {}

      reverse_iterator(int i, size_t n, T* data):m_i(n-1-i),m_n(n),m_data(data) {}

      void operator++() { m_i = ++m_i%m_n; }

      reverse_iterator operator+(int i) const { return reverse_iterator((m_i+i)%m_n,m_n,m_data); }

      void operator+=(int i) { m_i=(m_i+i)%m_n; }

      T& operator*() { return m_data[m_n-1-m_i]; }

      T* operator->() { return &(operator*()); }

      bool operator!=(const reverse_iterator& it) const { return m_i!=it.m_i; }

    private:

        int m_i;
        size_t m_n;
        T* m_data;

    };

    class const_reverse_iterator {

    public:

      const_reverse_iterator():m_i(0),m_n(0),m_data(nullptr) {}

      const_reverse_iterator(int i, size_t n, const T* data):m_i(n-1-i),m_n(n),m_data(data) {}

      void operator++() { m_i = ++m_i%m_n; }

      const_reverse_iterator operator+(int i) const { return const_reverse_iterator((m_i+i)%m_n,m_n,m_data); }

      void operator+=(int i) { m_i=(m_i+i)%m_n; }

      const T& operator*() const { return m_data[m_n-1-m_i]; }

      const T* operator->() const { return &(operator*()); }

      bool operator!=(const const_reverse_iterator& it) const { return m_i!=it.m_i; }

    private:

        int m_i;
        size_t m_n;
        const T* m_data;

    };

    //! Constructor.
    CRingBuffer(): m_data(),m_cursor(0) {}

    //! Constructor.
    explicit CRingBuffer(size_t n): m_data(n), m_cursor(0) {}

    //! Size.
    size_t size() const { return m_data.size(); }

    //! Pushback.
    void push_back(const T& x) { m_data[m_cursor] = x; m_cursor = (++m_cursor)%m_data.size(); }

    /*! \brief Back.
     *
     * This is the last element pushed into the container.
     *
     */
    const T& back() const { return m_data[(m_cursor+m_data.size()-1)%m_data.size()]; }

    //! Back.
    T& back() { return m_data[(m_cursor+m_data.size()-1)%m_data.size()]; }

    /*! \brief Front.
     *
     * Given a buffer length of \f$n\f$, the front is the element which was pushed \f$n-1\f$ times
     * before the back.
     *
     */
    const T& front() const { return m_data[m_cursor%m_data.size()]; }

    /*! \brief Front.
     *
     * Given a buffer length of \f$n\f$, the front is the element which was pushed \f$n-1\f$ times
     * before the back.
     *
     */
     T& front() { return m_data[m_cursor%m_data.size()]; }

    /*! \brief Creates forward iterator pointing to the "beginning".
     *
     * Beginning is always one element past the cursor.
     *
     */
    iterator begin() { return iterator((m_cursor)%m_data.size(),m_data.size(),&m_data[0]); }
    const_iterator begin() const { return const_iterator((m_cursor)%m_data.size(),m_data.size(),&m_data[0]); }

    /*! \brief Creates forward iterator pointing to the "end".
     *
     * TODO: This is suboptimal because it effectively shortens the buffer by
     * one element.
     *
     */
    iterator end() { return iterator((m_cursor+m_data.size()-1)%m_data.size(),m_data.size(),&m_data[0]); }
    const_iterator end() const { return const_iterator((m_cursor+m_data.size()-1)%m_data.size(),m_data.size(),&m_data[0]);  }

    //! Backward iterator.
    reverse_iterator rbegin() { return reverse_iterator((m_cursor+m_data.size()-1)%m_data.size(),m_data.size(),&m_data[0]); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator((m_cursor+m_data.size()-1)%m_data.size(),m_data.size(),&m_data[0]); }

    //! Backward iterator.
    reverse_iterator rend() { return reverse_iterator((m_cursor)%m_data.size(),m_data.size(),&m_data[0]); }
    const_reverse_iterator rend() const { return const_reverse_iterator((m_cursor)%m_data.size(),m_data.size(),&m_data[0]); }

    //! Reference to data.
    const std::vector<T>& data() { return m_data; }

    //! Access to cursor.
    const int& cursor() { return m_cursor; }

    //! Resize.
    void resize(size_t n) { m_data.resize(n); }

private:

    std::vector<T> m_data;
    int m_cursor;

};


}

#endif // RBUFFER_H
