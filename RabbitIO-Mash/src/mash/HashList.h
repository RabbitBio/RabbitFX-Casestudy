// Copyright © 2015, Battelle National Biodefense Institute (BNBI);
// all rights reserved. Authored by: Brian Ondov, Todd Treangen,
// Sergey Koren, and Adam Phillippy
//
// See the LICENSE.txt file included with this software for license information.

#ifndef HashList_h
#define HashList_h

#include "hash.h"
#include <vector>
#include <algorithm>

class HashList
{
public:
    
    HashList() {use64 = true;}
    HashList(bool use64new) {use64 = use64new;}
    
    hash_u at(int index) const;
    void clear();
    void resize(int size);
    void set32(int index, uint32_t value);
    void set64(int index, uint64_t value);
    void setUse64(bool use64New) {use64 = use64New;}
    int size() const {return use64 ? hashes64.size() : hashes32.size();}
    void sort();
    void push_back32(hash32_t hash) {hashes32.push_back(hash);}
    void push_back64(hash64_t hash) {hashes64.push_back(hash);}
    bool get64() const {return use64;}
	void merge(HashList & that)
	{
		if(use64){
    		std::vector<hash64_t> out64;
			std::merge(hashes64.begin(), hashes64.end(), 
					  that.hashes64.begin(), that.hashes64.end(),
                      std::back_inserter(out64));
			hashes64 = out64;
		}else{
    		std::vector<hash32_t> out32;
			std::merge(hashes32.begin(), hashes32.end(), 
					  that.hashes32.begin(), that.hashes32.end(),
                      std::back_inserter(out32));
			hashes32 = out32;

		}

		return;	
	}

public:
    
    bool use64;
    std::vector<hash32_t> hashes32;
    std::vector<hash64_t> hashes64;
};

#endif
