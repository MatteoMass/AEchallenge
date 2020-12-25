#pragma once

#include <algorithm>
#include <math.h>
#include <vector>
#include <stdint.h>
#include <iostream>
#include <sdsl/bit_vectors.hpp>

class MyIndex {
  size_t n;                             ///< The number of elements this index was built on.
  const std::vector<uint64_t> &data;    ///< The data vector this index was built on.
  int sizeBlock;
  std::vector<int> j;


public:

  MyIndex(const std::vector<uint64_t> &data) : n(data.size()), data(data) {
    int nBlock = pow(2,22);
    sizeBlock = ceil(double(data.back())/double(nBlock));
    j.resize(nBlock+1,-1);
    j.back() = data.size()-1;

    for(int x = 0; x < data.size(); ++x){
        int m = (data[x]/sizeBlock);
        if (j[m] != -1){
            ;
        }
        else{
            j[m] = x;
        }
    }
    
    int curr = 0;
    for(int x = 0; x < j.size(); ++x){
        if(j[j.size()-x] == -1){
            j[j.size()-x] = curr;
        }else{
            curr = j[j.size()-x];
        }
    }
  }
  

  /** Returns the element greater than or equal to a given key. */
  uint64_t nextGEQ(uint64_t key) const { // This function is required for the challenge
    int search = key/sizeBlock;
    int start = j[search];
    int end = j[search+1];
    return *std::lower_bound(data.begin()+start, data.begin()+end, key);

  }

  /** Returns the size of the index in bytes. */
  size_t size_in_bytes() const { // This function is required for the challenge
    return sizeof(this) + sizeof(int) * j.size();
    
  }
};