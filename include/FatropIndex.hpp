/**
 * @file FatropIndex.hpp
 * @author your name (you@domain.com)
 * @brief  index counter class
 * @version 0.1
 * @date 2021-11-09
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef FATROP_INDEX_INCLUDED
#define FATROP_INDEX_INCLUDED
using namespace std;
namespace fatrop
{
    class fatrop_index_counter
    {
    public:
        int next_index(const int size)
        {
            const int res = curr_index;
            curr_index += size;
            return res;
        }
        int get_size() const
        {
            return curr_index;
        }
    private:
        int curr_index = 0;
    };
}

#endif //FATROP_INDEX_INCLUDED