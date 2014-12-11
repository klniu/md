#ifndef SYMMETRY_H
#define SYMMETRY_H
#include <openbabel/mol.h>
#include <openbabel/isomorphism.h>
#include <openbabel/query.h>
#include <set>

using namespace OpenBabel;
typedef std::pair<unsigned, unsigned> upair;
typedef std::vector<std::pair<unsigned, unsigned> > vupair;
typedef std::vector<std::set<unsigned> > vuset;
typedef std::vector<unsigned> uvector;
typedef std::set<unsigned> uset;
/**
 * @brief Get the symmetrys of a molecule. In the symmetrys, all the symmetry atom will be in one vector.
 *
 */
class Symmetry
{
public:
/**
 * @brief
 *
 * @param mol The molecule.
 * @param timeout The timeout of processing.
 */
    Symmetry(OBMol *mol, unsigned timeout);
    ~Symmetry();

    /**
     * @brief Get the symmetrys of the molecule.
     *
     * @return vuset
     */
    vuset &getSymmetry();

private:
    /**
     * @brief Traversal all the item in \a amlist, finding the tree from \a items. e.g items = [1,2], amlist = [[1,3], [2,4], [4,5], [6,7]], the \a sameSymmetry will be [1,2,3,4,5]. Like a tree.
     *
     * @param items A unsigned std::Set Include the items which are to be found.
     * @param amlist Where to find the item
     * @param sameSymmetry The result.
     */
    void findSameSymmetry(uset items, vupair &amlist, vupair &sameSymmetry);

private:
    OBIsomorphismMapper *mapper;
    unsigned timeout;
    vupair symmetryPairs; /**< TODO The initial symmetry pairs geting from mapperfunctor*/
    vuset symmetrys; /**< TODO The last result*/
};

#endif // SYMMETRY_H
