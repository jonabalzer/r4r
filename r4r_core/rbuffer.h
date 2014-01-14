#ifndef R4RRBUFFER_H
#define R4RRBUFFER_H

#include <stdlib.h>
#include <deque>



namespace R4R {

/*! a simple ring buffer implementation
 *
 * This class inherits from std::deque so that the contents of the container can be
 * accessed through STL iterators of the base class. One should also inhibit any
 * operations that change the size of the ring buffer and that are available through
 * the implementation of the base class.
 *
 */
template<class T>
class CRingBuffer: public std::deque<T> {

public:

    //! Delete standard constructor to make sure the size of the ring buffer is specified.
    CRingBuffer() = delete;

    //! Constructor.
    explicit CRingBuffer(size_t n): std::deque<T>(n) {}

    //! Insertion operation.
    void push(const T& x) {

        std::deque<T>::pop_back();
        std::deque<T>::push_front(x);

    }

private:

};


}

#endif // RBUFFER_H
