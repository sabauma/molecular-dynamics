/**
 * This class provides the functionality of a three dimensional cell list. It
 * provides the ability to quickly find neighboring particles in a three
 * dimensional system. This is accomplished using a three dimensional
 * array to store the contents of each voxel in the system. In that manner,
 * neighbors can be quickly found by examining the current and adjacent
 * voxels.
 *
 * @author Spenser Bauman
 */

#ifndef __CELL_LIST_CELL_LIST_FAST_H__
#define __CELL_LIST_CELL_LIST_FAST_H__

#include <boost/multi_array.hpp>

#include "../Vector/Vector.hpp"

class CellList
{
private:

    typedef boost::multi_array<int, 3> CellData;
    typedef boost::multi_array<int, 1> LinkedList;
    typedef boost::multi_array<IntVector, 1> CellIndex;

    int Particles;          // The number of particles in the system
    IntVector Size;

    CellData Cells;         // 3D array of cell contents.
    CellIndex Index;        // An index of which cell each particle is in.
    LinkedList Neighbors;   // A list of the neighbors in each cell.

    /**
     * Traverses the list beginning at index <code>start</code> in the
     * Neighbors array. This advances along the list until the element just
     * before the element <code>ps</code> and returns that index. This
     * makes. This is mostly used for deleting elements from
     * the list.
     *
     *
     * @param start The head of the list to being searching.
     * @param ps The particle to find in the list.
     * @return The index of the element in the list before the particle
     *         ps. -1 is returned if start == ps (i.e. ps is the first
     *         element of the list).
     */
    int FindListPosition(int start, const int ps);

public:

    class NeighborCellIterator
    {
    private:

        IntVector LowBounds;
        IntVector HighBounds;

        CellList& Parent;                          // Reference to the current cell list

        const IntVector Cell;                      // Cell to find neighbors for
        IntVector CurrentCell;
        int CurrentParticle;

    public:

        // Constructor for the inner cell
        NeighborCellIterator(CellList& parent, const IntVector& cell,
                             const IntVector& low, const IntVector& high);

        // De-reference Operation
        int operator*() const;

        // Have we exhaused the iterator?
        bool HasNext() const;

        // Increment the iterator
        NeighborCellIterator& operator++();
    };

    /**
     * Constructor
     */
    CellList(const int ps, const int xs, const int ys, const int zs);

    /**
     * Get the particle's position in the grid space
     */
    const IntVector& ParticlePosition(const int ps) const;

    /**
     * Return the particle position through a referenced IntVector.
     */
    void ParticlePosition(const int ps, IntVector& v) const;

    /**
     * Moves a particle the set number of cells in each direction.
     *
     * @param ps The index of the particle to move.
     * @param mv The displacement in each direction.
     */
    void MoveParticle(const int ps, const IntVector& mv);

    /**
     * Sets the position of the current particle.
     *
     * @param particle The index of the particle being set.
     * @param x The x coordinate of the particle.
     * @param y The y coordinate of the particle.
     * @param z The z coordinate of the particle.
     */
    void SetParticle(const int ps, const int x, const int y, const int z);


    /**
     * Sets the position of the current particle. Functions as a convenient
     * adaptor as we often wish to position the particles based on a pre
     * computed vector.
     *
     * @param particle The index of the particle being set.
     * @param v Vector representing the coordinates of the particle.
     */
    void SetParticle(const int particle, const IntVector& v);

    /**
     * Create an iterator into the neighboring cells of a given particle.
     * Will include the current particle.
     *
     * @param particle The particle to get the neighbors for.
     * @param low The low bounds along each axis for traversal. This is
     *            relative to the position of the particle with the given
     *            index.
     * @param high The high bound along each axis.
     */
    NeighborCellIterator BeginNeighborCell(const int particle,
                                           const IntVector& low,
                                           const IntVector& high);

    /**
     * Resize the cell list to a new volume. This will invalidate and clear
     * out the contents of the <code>CellList</code>. Any iterators into this
     * object will be invalidated once this operaton is performed.
     *
     * @param ps The number of particles in the new system.
     * @param xs The number of partitions along the x-axis.
     * @param ys The number of partitions along the y-axis.
     * @param zs The number of partitions along the z-axis.
     */
    void Resize(const int ps, const IntVector& ex);

    /**
     * Clear out the contents of the <code>CellList</code>. This will
     * invalidate all of the particle placements.
     */
    void Clear();

    /**
     * Get the number of partitions along each dimension.
     *
     * @param i The index of the dimension (0 = x, 1 = y, 2 = z).
     * @return The number of partitions along the given dimension.
     */
    inline unsigned long Dims(const int i) const
    {
        assert(i >= 0);
        assert((unsigned) i < Cells.num_dimensions());
        return Cells.shape()[i];
    }
};

#endif  /* __CELL_LIST_H__ */
