#include "mapperfunctor.h"

bool MapperFunctor::operator()(OBIsomorphismMapper::Mapping &map)
{
    std::pair<unsigned, unsigned> p2;
    for (std::vector<std::pair<unsigned, unsigned> >::iterator p = map.begin(); p != map.end(); ++p){
        if (p->first == p->second){
            continue;
        }else{
            // sort
            if (p->first > p->second)
                p2 = std::make_pair(p->second + 1, p->first + 1);
            else
                p2 = std::make_pair(p->first + 1, p->second + 1);

            std::vector<std::pair<unsigned, unsigned> >::iterator exist = std::find(sym_data.begin(), sym_data.end(), p2);

            if (exist == sym_data.end())
                sym_data.push_back(p2);
        }
        //for (std::vector<std::pair<unsigned, unsigned> >::iterator i = sym_data.begin(); i != sym_data.end(); ++i){
        //    std::cout << i->first << i->second << std::endl;
        //}
    }
    // code to handle found isomorphism...
    // examples: - store the mapping
    //           - filter mappings
    //           - use the found map in some way


    // continue mapping
    return false;
}
