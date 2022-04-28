#pragma once

#include <map>
#include <vector>

#include "parallel/OpenMP.h"

namespace templatetools
{
    template <typename Key, typename Value>
    void InsertIntoMap(Key& key, Value& value, std::map<Key,Value>& map) 
    {
        typename std::map<Key,Value>::iterator it = map.find(key);

        if (it == map.end())
        {
            map.insert(std::make_pair(key, value));
        }

        return;
    }

    // maps from something to a vector of something
    template <typename Key, typename Value>
    void InsertIntoVectorMap(Key& k, Value& v, std::map<Key,std::vector<Value>>& map) 
    {
        typename std::map<Key,std::vector<Value>>::iterator it = map.find(k);

        if (it == map.end())
        {
            std::vector<Value> tempVec = {v};
            map.insert(std::make_pair(k,tempVec));
        }
        else
        {
            it -> second.push_back(v);
        }

        return;
    }
}