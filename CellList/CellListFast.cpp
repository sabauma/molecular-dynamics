
#include <algorithm>
#include <assert.h>
#include <boost/multi_array.hpp>
#include <stdio.h>
#include <iostream>

#include "CellListFast.h"
#include "../Vector/Vector.hpp"

/**
 * Incrmenets the cell position by one assuming lower and upper bounds
 * on each axis for a carry over effect.
 */
static inline void
increment(IntVector& currentCell, const IntVector& low, const IntVector high)
{
    const int end = (int) currentCell.size() - 1;
    // Increment lowest value
    ++currentCell[end];

    // Perform carry over
    for (int i = end; i > 0 && currentCell[i] > high[i]; --i)
    {
        currentCell[i] = low[i];
        ++currentCell[i-1];
    }
}

int
CellList::FindListPosition(int start, const int ps)
{
    assert(start >= 0);
    assert(start < (int) this->Neighbors.size() - 1);

    if (start == ps) return -1;

    // Traverse the list until we find the one before the index
    // we wish to remove
    while (this->Neighbors[start] != ps)
    {
        assert(this->Neighbors[start] != -1);
        start = this->Neighbors[start];
    }

    return start;
}

/**
 * Fixes a 3d vector such that it "wraps" around a space with dimensions
 * specified by x, y, and z.
 *
 * @param v The vector to be fixed.
 * @param x The size of the x dimension.
 * @param y The size of the y dimension.
 * @param z The size of the z dimension.
 */
CellList::CellList(const int ps, const int xs, const int ys, const int zs) :
    Particles(ps), Size(xs, ys, zs),
    Cells(boost::extents[xs][ys][zs]),
    Index(boost::extents[ps]),
    Neighbors(boost::extents[ps])
{
    IntVector* index_start = Index.data();
    IntVector* index_end   = index_start + Index.num_elements();

    int* neighbors_start   = this->Neighbors.data();
    int* neighbors_end     = neighbors_start + this->Neighbors.num_elements();

    int* cells_start       = this->Cells.data();
    int* cells_end         = cells_start + this->Cells.num_elements();

    // Fill with default values
    std::fill(index_start, index_end, -1);
    std::fill(neighbors_start, neighbors_end, -1);
    std::fill(cells_start, cells_end, -1);
}

const IntVector&
CellList::ParticlePosition(const int ps) const
{
    return Index[ps];
}

void
CellList::ParticlePosition(const int ps, IntVector& v) const
{
    v = Index[ps];
}

void
CellList::MoveParticle(const int ps, const IntVector& mv)
{
    assert(ps >= 0 && ps < Particles);

    // The cell position of this particle.
    IntVector next_place(Index[ps] + mv);
    this->SetParticle(ps, next_place);
}

/**
 * Get an iterator into the neighboring cells of a particle, including the
 * one the current particle is occupying.
 *
 * @param particle The index of the particle to search around.
 */
CellList::NeighborCellIterator
CellList::BeginNeighborCell(
        const int particle,
        const IntVector& low,
        const IntVector& high)
{
    // Must operate on a valid particle index.
    assert(particle >= 0 && particle < Particles);
    // Get its cell
    IntVector& cell = this->Index[particle];
    // Cannot get particles around an uninitialized particle.
    assert(all_elements(cell != -1));

    return CellList::NeighborCellIterator(*this, cell, low, high);
}

void
CellList::SetParticle(const int particle, const int x, const int y, const int z)
{
    const IntVector v(x, y, z);
    this->SetParticle(particle, v);
}

void
CellList::SetParticle(const int particle, const IntVector& v)
{
    assert(particle >= 0 && particle < Particles);
    assert(all_elements(v >= 0));
    assert(all_elements(v < Size));

    IntVector& cell = this->Index[particle];

    assert(all_elements(cell == -1) || all_elements(cell != -1));

    // The particle must be in a cell, if it is not the initial value.
    // Therefore, delete it from the current list it is in.
    if (cell[0] != -1 && cell[1] != -1 && cell[2] != -1)
    {
        int& bucket = this->Cells[cell[0]][cell[1]][cell[2]];

        // Particle is in the system
        const int c = this->FindListPosition(bucket, particle);

        // If the element is not the first one on the list, then we need to
        // fix up the list. Clear out the list if this is the first element.
        if (c != -1)
        {
            this->Neighbors[c] = this->Neighbors[c] != -1
                               ? this->Neighbors[this->Neighbors[c]]
                               : -1;
        }
        else
        {
            // The next item becomes the new particle.
            bucket = this->Neighbors[particle];
            // Invalidate this particle.
            this->Neighbors[particle] = -1;
        }
    }

    cell = v;

    // Place the particle in the list corresponding to the cell it moved into.
    // Particle gets placed at the head of the list.
    int& ref = Cells[cell[0]][cell[1]][cell[2]];
    int tmp = ref;
    ref = particle;
    this->Neighbors[particle] = tmp;
}

/**
 * Create an iterator into the cells adjacent to the current one (the current
 * one is treated as trivially adjacent)
 *
 * @param parent The cell list that the iterator is being created for.
 * @param cell The central cell we are searching from.
 */
CellList::NeighborCellIterator::NeighborCellIterator(
        CellList& parent, const IntVector& cell,
        const IntVector& low, const IntVector& high) :
    Parent(parent), Cell(cell)
{
    // Determine the bounds of space to explore.
    for (int i = 0; i < 3; ++i)
    {
        CurrentCell[i] = LowBounds[i] = std::max(-Cell[i], low[i]);
        HighBounds[i]  = std::min(Parent.Size[i] - Cell[i] - 1, high[i]);
    }

    IntVector next(Cell + LowBounds);
    assert(all_elements(next >= 0));

    // Edge case where everything is out of bounds. This means that the
    // initial low bounds are out of range, so we need not continue searching.
    if (any_elements(next >= Parent.Size))
    {
        this->CurrentParticle = -1;
        return;
    }
    else
    {
        this->CurrentParticle = Parent.Cells[next[0]][next[1]][next[2]];
    }

    // Advance to the first non-empty cell.
    while (this->CurrentParticle == -1)
    {
        increment(this->CurrentCell, this->LowBounds, this->HighBounds);

        if (CurrentCell[0] > HighBounds[0])
        {
            break;
        }

        next = this->Cell + this->CurrentCell;

        this->CurrentParticle = Parent.Cells[next[0]][next[1]][next[2]];
    }
}

/**
 * Advances the iterator one position.
 * @return a reference to the current iterator.
 */
CellList::NeighborCellIterator&
CellList::NeighborCellIterator::operator++()
{
    assert(this->HasNext());
    assert(this->CurrentParticle != Parent.Neighbors[this->CurrentParticle]);
    IntVector next;

    // Advance to the next particle in the current list.
    this->CurrentParticle = Parent.Neighbors[this->CurrentParticle];

    // If we are out of particles in the current cell, search
    // through adjacent ones until we find another populated cell.
    while (this->CurrentParticle == -1)
    {
        increment(this->CurrentCell, this->LowBounds, this->HighBounds);

        if (CurrentCell[0] > HighBounds[0])
        {
            break;
        }

        assert(this->CurrentParticle >= -1 &&
               this->CurrentParticle < Parent.Particles);

        next = this->Cell + this->CurrentCell;

        this->CurrentParticle = Parent.Cells[next[0]][next[1]][next[2]];
    }


    return *this;
}

/**
 * @return The current value the iterator is focused on.
 */
int
CellList::NeighborCellIterator::operator*() const
{
    return this->CurrentParticle;
}

/**
 * @return Returns true if the iterator can be safely dereferenced.
 */
bool
CellList::NeighborCellIterator::HasNext() const
{
    // The iterator is valid if the current adjacent cell is within bounds
    // and the iterator is not exhausted.
    return this->CurrentParticle != -1;
}

void
CellList::Resize(const int ps, const IntVector& ex)
{
    Cells.resize(boost::extents[0][0][0]);
    Index.resize(boost::extents[0]);
    Neighbors.resize(boost::extents[0]);

    // Resize to the desired size.
    Cells.resize(boost::extents[ex[0]][ex[1]][ex[2]]);
    Index.resize(boost::extents[ps]);
    Neighbors.resize(boost::extents[ps]);

    IntVector* index_start = Index.data();
    IntVector* index_end   = index_start + Index.num_elements();

    int* neighbors_start   = this->Neighbors.data();
    int* neighbors_end     = neighbors_start + this->Neighbors.num_elements();

    int* cells_start       = this->Cells.data();
    int* cells_end         = cells_start + this->Cells.num_elements();

    // Fill with default values
    std::fill(index_start, index_end, -1);
    std::fill(neighbors_start, neighbors_end, -1);
    std::fill(cells_start, cells_end, -1);

    this->Particles = ps;
    this->Size      = ex;
}

void CellList::Clear()
{
    IntVector* index_start = Index.data();
    IntVector* index_end   = index_start + Index.num_elements();

    int* neighbors_start   = this->Neighbors.data();
    int* neighbors_end     = neighbors_start + this->Neighbors.num_elements();

    int* cells_start       = this->Cells.data();
    int* cells_end         = cells_start + this->Cells.num_elements();

    // Fill with default values
    std::fill(index_start, index_end, -1);
    std::fill(neighbors_start, neighbors_end, -1);
    std::fill(cells_start, cells_end, -1);
}
