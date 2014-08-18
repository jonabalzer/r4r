//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, Jonathan Balzer
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
//////////////////////////////////////////////////////////////////////////////////

#ifndef R4RUNIONFIND_H
#define R4RUNIONFIND_H

#include <vector>
#include <list>
#include <iostream>
#include <stdlib.h>

namespace R4R {

template<class T>
class CUnionFindContainer {

public:

    /*! \brief list that represents each element of the partitioon
     *
     */
    class CUFList {

    friend class CUnionFindContainer<T>;

    public:

        /*!
         * \brief list node
         */
        class CUFListNode {

        friend class CUFList;

        public:

            CUFListNode(const T& data):m_head(nullptr),m_next(nullptr),m_data(data) {}

            CUFListNode():m_head(nullptr),m_next(nullptr),m_data() {}

            void PrintChildren(std::ostream& out) const {
                out << m_data << std::endl;
                if(m_next!=nullptr)
                    m_next->PrintChildren(out);
            }

            void SetData(const T& data) { m_data = data; }
            T& GetData() const { return m_data; }


        private:
            CUFListNode* m_head;
            CUFListNode* m_next;
            T m_data;

        };


        //! Constructor
        CUFList():m_size(0),m_head(nullptr),m_tail(nullptr) {}

        CUFList(CUFListNode* node):m_size(0),m_head(nullptr),m_tail(nullptr) { this->PushBack(node);}


        /*! \brief Push a value into the list.
         *
         * Do node allocation outside of class.
         *
         *
         */
        void PushBack(CUFListNode* node) {

            // is this the first node
            if(m_size==0) {

                node->m_head = node;
                m_head = node;
                m_tail = node;

            }
            else {
                node->m_head = m_head;
                m_tail->m_next = node;
                m_tail = node;
            }

            // increment cursor
            m_size++;

        }

        size_t Size() const { return m_size; }

        void Print(std::ostream& out) const {

            if(m_head!=nullptr)
                m_head->PrintChildren(out);

        }

    private:

        size_t m_size;
        CUFListNode* m_head;
        CUFListNode* m_tail;

    };

    typedef typename CUFList::CUFListNode CUFListNode;

    //! Constructor.
    CUnionFindContainer():m_partitions() {}

    //! Returns the represenative of a list.
    const CUFListNode* Find(const CUFListNode* x) { return x->m_head; }

    /*! \brief Creates a new set.
     *
     */
    void MakeSet(CUFListNode* node) {

        // only add if it is not in any set
        if(node->m_head!=nullptr) {
            CUFList newlist(node);
            m_partitions.push_back(newlist);
        }




    }

    //! Union operation.
    void Union() {}




private:

    std::vector<CUFList> m_partitions;

};

//template<class T>


//    CUnionFindContainer():m_partitions() {}

//    void MakeSet(const T& x) {

//        CUFList<T> partition;
//        partition.push_back(x);
//        m_partitions.push_back(partition);

//    }



//private:

//    std::vector<CUFList<T> > m_partitions;

//};

} // end of namespace

#endif // UNIONFIND_H
