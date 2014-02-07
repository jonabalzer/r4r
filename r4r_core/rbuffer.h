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

    //! Delete standard constructor to make sure the size of the ring buffer is specified.
    CRingBuffer() = delete;
    CRingBuffer(const CRingBuffer<T>& x) = delete;
    CRingBuffer<T>& operator=(const CRingBuffer<T>& x) = delete;

    //! Constructor.
    explicit CRingBuffer(size_t n): m_data(n), m_cursor(0) {}

    //! Size.
    size_t size() const { return m_data.size(); }

    //! Pushback.
    void push_back(const T& x) { m_data[m_cursor] = x; m_cursor = (++m_cursor)%m_data.size(); }

    //! Back
    const T& back() const { return m_data[(m_cursor+m_data.size()-1)%m_data.size()]; }

    //! Back
    T& back() { return m_data[(m_cursor+m_data.size()-1)%m_data.size()]; }

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
