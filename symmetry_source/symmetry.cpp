#include "symmetry.h"
#include "mapperfunctor.h"

Symmetry::Symmetry(OBMol *mol, unsigned timeout)
{

//    OBGraphSym graphsym(mol);
//    unsigned a = 0;
//    std::vector<unsigned> symmetry_classes(a);
//    graphsym.GetSymmetry(symmetry_classes);
//    std::vector<OBIsomorphismMapper::Mapping> automorphisms;
//    MapperFunctor mapperFunctor(automorphisms);
//    FindAutomorphisms(mapperFunctor, mol, symmetry_classes);


    std::vector<OBIsomorphismMapper::Mapping> automorphisms;
    MapperFunctor mapperFunctor(automorphisms);
    OBQuery *query = CompileMoleculeQuery(mol);
    mapper = OBIsomorphismMapper::GetInstance(query);
    mapper->SetTimeout(timeout);
    mapper->MapGeneric(mapperFunctor, mol);
    symmetryPairs = mapperFunctor.getSymmetry();
}

Symmetry::~Symmetry()
{
    delete mapper;
}

vuset & Symmetry::getSymmetry()
{
    vupair pairs(symmetryPairs);
    while (pairs.size() > 0){
        uset items;
        items.insert(pairs.at(0).first);
        items.insert(pairs.at(0).second);
        vupair sameSymmetryPairs;
        findSameSymmetry(items, symmetryPairs, sameSymmetryPairs);

        uset foundItems;
        for (vupair::iterator j = sameSymmetryPairs.begin(); j != sameSymmetryPairs.end(); ++j){
            foundItems.insert(j->first);
            foundItems.insert(j->second);
            vupair::iterator findit = std::find(pairs.begin(), pairs.end(), *j);
            if (findit != pairs.end()){
                pairs.erase(findit);
            }
        }
        symmetrys.push_back(foundItems);
    }
    return symmetrys;
}

void Symmetry::findSameSymmetry(uset items, vupair &amlist, vupair &sameSymmetry)
{
    uset usedItems;
    while (items.size() > 0){
        unsigned usingItem = *(items.begin());
        usedItems.insert(usingItem);
        for (vupair::iterator i = amlist.begin(); i != amlist.end(); ++i){
            if (usingItem == i->first || usingItem == i->second){
                sameSymmetry.push_back(*i);
                items.insert(i->first);
                items.insert(i->second);
            }
        }
        for (uset::iterator i = usedItems.begin(); i != usedItems.end(); ++i){
            items.erase(*i);
        }
    }
}
