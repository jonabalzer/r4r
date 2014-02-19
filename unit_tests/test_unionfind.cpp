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

#include <gtest/gtest.h>
#include <string>
#include "unionfind.h"

using namespace R4R;
using namespace std;

class CUnionFindContainerTests: public testing::Test {

protected:

    CUnionFindContainerTests():
        m_free_store(10),
        m_list()
        {}

    virtual void SetUp() {

        for(size_t i=0; i<m_free_store.size(); i++) {
            // allocate here!!!
            //m_free_store[0].SetData(0.1);
            //m_free_store[1].SetData(0.2);
            //m_free_store[2].SetData(0.3);
            //m_free_store[3].SetData(0.4);
            //m_list.PushBack(&m_free_store[0]);
            // put it into a list
        }
    }

    virtual void TearDown() {



    }


    std::vector<CUnionFindContainer<float>::CUFListNode*> m_free_store;
    CUnionFindContainer<float>::CUFList m_list;                             //!< list used in union-find data structure
    //CUnionFindContainer<float> m_container;

};

TEST_F(CUnionFindContainerTests, ListTests) {

    // checking addition and iteration
//    stringstream ss, sse;
//    m_list.Print(ss);
//    sse << 0.1 << endl;
//    sse << 0.2 << endl;
//    sse << 0.3 << endl;
//    sse << 0.4 << endl;
//    EXPECT_EQ(ss.str(),sse.str());

    // add seven more items and check
//    m_list.PushBack(0.5);
//    m_list.PushBack(0.6);
//    m_list.PushBack(0.7);
//    m_list.PushBack(0.8);
//    m_list.PushBack(0.9);
//    m_list.PushBack(1.0);
//    m_list.PushBack(1.1);
//    EXPECT_EQ(11,m_list.Size());

}

