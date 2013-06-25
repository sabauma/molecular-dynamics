/**
 *
 * @author Spenser Bauman
 */

#ifndef __CELL_LIST_H__
#define __CELL_LIST_H__

#include <boost/multi_array.hpp>
#include <set>

#include "../Vector/Vector.hpp"

class CellList
{
private:

    typedef boost::multi_array<std::set<int>, 3> CellData;
    typedef boost::multi_array<IntVector, 1> CellIndex;

    int Particles;    // The number of particles in the system
    IntVector Size;

    CellData Cells;         // 3D array of cell contents.
    CellIndex Index;        // An index of which cell each particle is in.

public:

    // Get an iterator of all adjacent particles.
    class CurrentCellIterator
    {
    private:

        const IntVector Cell;                 // The cell the iterator is in.
        std::set<int>::const_iterator Iter;   // Iterator into the cell
        std::set<int>::const_iterator End;    // End iterator for the cell

    public:
        // Constructor for the inner cell
        CurrentCellIterator(CellList& parent,
                            const IntVector cell);

        // De-reference Operation
        int operator*() const;

        // Does the iterator have more to give?
        bool HasNext() const;

        // Increment the iterator
        CurrentCellIterator& operator++();
    };

    class NeighborCellIterator
    {
    private:

        IntVector LowBounds;
        IntVector HighBounds;

        std::set<int>::const_iterator CurrentIter; // Iterator into the current cell
        std::set<int>::const_iterator EndIter;     // End iterator for the current cell
        CellList& Parent;                    // Reference to the current cell list

        const IntVector Cell;                // Cell to find neighbors for
        IntVector CurrentCell;

    public:
        // Constructor for the inner cell
        NeighborCellIterator(CellList& parent, const IntVector& cell);

        // De-reference Operation
        int operator*();

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
     * Get an iterator into the same cell as the current particle.
     * Will include the current particle.
     *
     * @param particle The index of the particle to search around.
     */
    CurrentCellIterator BeginSameCell(const int particle);

    /**
     * Create an iterator into the neighboring cells of a given particle.
     * Will include the current particle.
     *
     * @param particle The particle to get the neighbors for.
     */
    NeighborCellIterator BeginNeighborCell(const int particle);

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
    void Resize(const int ps, const int xs, const int ys, const int zs);

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
