#include <iostream>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include "symmetry.h"

using namespace OpenBabel;

typedef std::pair<unsigned, unsigned> upair;
typedef std::vector<std::pair<unsigned, unsigned> > vupair;
typedef std::vector<std::set<unsigned> > vuset;
typedef std::vector<unsigned> uvector;
typedef std::set<unsigned> uset;
int main(int argc,char **argv)
{
    if(argc<3)
    {
        std::cout << "Usage: Symmetry MoleculeFileFormat InputFileName\n";
        return 1;
    }

    std::ifstream ifs(argv[2]);
    if(!ifs)
    {
        std::cout << "Cannot open input file\n";
        return 1;
    }
    OBConversion obconversion;
    OBMol mol;
    if (!obconversion.SetInFormat(argv[1])){
        std::cout << "Formats not available\n";
        return 1;
    }
    if (obconversion.ReadFile(&mol, argv[2])){
        // OBGraphSym graphsym(&mol);
        // unsigned a = 0;
        // std::vector<unsigned> symmetry_classes(a);
        // graphsym.GetSymmetry(symmetry_classes);
        // std::vector<OBIsomorphismMapper::Mapping> automorphisms;
        // SymmetryMapFuctor symmetryMapFuctor(automorphisms);
        // FindAutomorphisms(symmetryMapFuctor, &mol, symmetry_classes);

        Symmetry symmetry(&mol, 36000);
        vuset symmetrys = symmetry.getSymmetry();
        std::string result = "[";
        for (vuset::iterator i = symmetrys.begin(); i != symmetrys.end(); ++i){
            result += std::to_string(*(i->begin()));
            for (uset::iterator j = ++(i->begin()); j != i->end(); ++j){
                result += ", " + std::to_string(*j);
            }
            result += "], [";
        }
        result.pop_back();
        result.pop_back();
        result.pop_back();
        result += "\n";
        result += "Number: " + std::to_string(symmetrys.size()) + "\n";
        std::cout << result;
    }

    return 0;
}
