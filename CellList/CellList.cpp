// TODO: Make the changes to the neighbor cell iterators to support ranges.
#include <assert.h>
#include <boost/multi_array.hpp>
#include <set>
#include <stdio.h>
#include <iostream>

#include "CellList.h"

static const int ADJACENT_CELLS = 27;

static IntVector ADJACENCY_LIST[ADJACENT_CELLS] =
{
    IntVector(-1 , -1 , -1 ) ,
    IntVector(-1 , -1 , 0  ) ,
    IntVector(-1 , -1 , 1  ) ,

    IntVector(-1 , 0  , -1 ) ,
    IntVector(-1 , 0  , 0  ) ,
    IntVector(-1 , 0  , 1  ) ,

    IntVector(-1 , 1  , -1 ) ,
    IntVector(-1 , 1  , 0  ) ,
    IntVector(-1 , 1  , 1  ) ,


    IntVector(0  , -1 , -1 ) ,
    IntVector(0  , -1 , 0  ) ,
    IntVector(0  , -1 , 1  ) ,

    IntVector(0  , 0  , -1 ) ,
    IntVector(0  , 0  , 0  ) ,
    IntVector(0  , 0  , 1  ) ,

    IntVector(0  , 1  , -1 ) ,
    IntVector(0  , 1  , 0  ) ,
    IntVector(0  , 1  , 1  ) ,


    IntVector(1  , -1 , -1 ) ,
    IntVector(1  , -1 , 0  ) ,
    IntVector(1  , -1 , 1  ) ,

    IntVector(1  , 0  , -1 ) ,
    IntVector(1  , 0  , 0  ) ,
    IntVector(1  , 0  , 1  ) ,

    IntVector(1  , 1  , -1 ) ,
    IntVector(1  , 1  , 0  ) ,
    IntVector(1  , 1  , 1  ) ,
};

CellList::CellList(const int ps, const int xs, const int ys, const int zs)
    : Particles(ps), Size(xs, ys, zs)
    , Cells(boost::extents[xs][ys][zs])
    , Index(boost::extents[ps])
{
    IntVector* start            = Index.data();
    IntVector const * const end = start + Index.num_elements();

    for (; start != end; ++start)
    {
        (*start)[0] = (*start)[1] = (*start)[2] = -1;
    }
}

const IntVector&
CellList::ParticlePosition(const int ps) const
{
    return Index[ps];
}

void
CellList::MoveParticle(const int ps, const IntVector& mv)
{
    assert(ps >= 0 && ps < Particles);

    // The cell position of this particle.
    IntVector& cur = Index[ps];

    // Remove it from the current cell.
    Cells[cur[0]][cur[1]][cur[2]].erase(ps);

    // Calculate the particles new position with a wrapped boundary.
    cur += mv;

    // Place the particle in its new cell.
    Cells[cur[0]][cur[1]][cur[2]].insert(ps);
}

CellList::CurrentCellIterator
CellList::BeginSameCell(const int particle)
{
    assert(particle >= 0 && particle < Particles);
    IntVector& cell = this->Index[particle];
    return CellList::CurrentCellIterator(*this, cell);
}

/**
 * Get an iterator into the neighboring cells of a particle, including the
 * one the current particle is occupying.
 *
 * @param particle The index of the particle to search around.
 */
CellList::NeighborCellIterator
CellList::BeginNeighborCell(const int particle)
{
    assert(particle >= 0 && particle < Particles);
    IntVector& cell = this->Index[particle];
    return CellList::NeighborCellIterator(*this, cell);
}

void
CellList::SetParticle(const int particle, const int x, const int y, const int z)
{
    assert(particle >= 0 && particle < Particles);
    assert(all_elements(IntVector(x, y, z) >= 0));
    assert(all_elements(IntVector(x, y, z) < Size));

    IntVector& cell = this->Index[particle];

    if (all_elements(cell == -1))
    {
        // Particle not in the system yet
        cell[0] = x;
        cell[1] = y;
        cell[2] = z;
    }
    else
    {
        // Particle is in the system
        this->Cells[cell[0]][cell[1]][cell[2]].erase(particle);
        cell[0] = x;
        cell[1] = y;
        cell[2] = z;
    }

    this->Cells[x][y][z].insert(particle);
}

void
CellList::SetParticle(const int particle, const IntVector& v)
{
    assert(v.size() == 3);
    this->SetParticle(particle, v[0], v[1], v[2]);
}

/**
 * Construct an iterator into the current cell.
 *
 * @param parent The cell list being iterated through.
 * @param cell The cell in the parent being iterated through.
 */
CellList::CurrentCellIterator::CurrentCellIterator(
        CellList& parent, const IntVector cell) : Cell(cell)
{
    std::set<int>& ref = parent.Cells[Cell[0]][Cell[1]][Cell[2]];
    Iter = ref.begin();
    End  = ref.end();
}

/**
 * Advance the iterator one position.
 *
 * @return A reference to the current iterator.
 */
CellList::CurrentCellIterator&
CellList::CurrentCellIterator::operator++()
{
    assert(Iter != End);
    ++Iter;
    return *this;
}

/**
 * Return the current focus of the iterator.
 */
int
CellList::CurrentCellIterator::operator*() const
{
    assert(this->HasNext());
    return *Iter;
}

/**
 * Checks that the iterator is safe to dereference.
 *
 * @return True if the iterator is safe to dereference.
 */
bool
CellList::CurrentCellIterator::HasNext() const
{
    return Iter != End;
}

/**
 * Create an iterator into the cells adjacent to the current one (the current
 * one is treated as trivially adjacent)
 *
 * @param parent The cell list that the iterator is being created for.
 * @param cell The central cell we are searching from.
 */
CellList::NeighborCellIterator::NeighborCellIterator(
        CellList& parent, const IntVector& cell)
    : Parent(parent), Cell(cell)
{
    for (int i = 0; i < 3; ++i)
    {
        CurrentCell[i] = LowBounds[i] = Cell[i] == 0 ? 0 : -1;
        HighBounds[i] = Cell[i] == Parent.Size[i] - 1 ? 0 : 1;
    }

    std::set<int>& ref = Parent.Cells[Cell[0] + LowBounds[0]]
                                     [Cell[1] + LowBounds[1]]
                                     [Cell[2] + LowBounds[2]];
    this->CurrentIter  = ref.begin();
    this->EndIter      = ref.end();

    // Advance to the first non-empty cell.
    while (CurrentIter == EndIter)
    {
        // Increment x
        ++CurrentCell[2];

        if (CurrentCell[2] > HighBounds[2]) // Carry over if x is out of bounds
        {
            ++CurrentCell[1];
            CurrentCell[2] = LowBounds[2];

            if (CurrentCell[1] > HighBounds[1]) // Carry over if y is out of bounds
            {
                ++CurrentCell[0];
                CurrentCell[1] = LowBounds[1];

                if (CurrentCell[0] > HighBounds[0]) // End if z is out of bounds.
                {
                    break;
                }
            }
        }

        const IntVector next(this->Cell + this->CurrentCell);

        std::set<int>& ref = Parent.Cells[next[0]][next[1]][next[2]];

        if (!ref.empty())   // Only get the iterator if they are not empty.
        {
            CurrentIter = ref.begin();
            EndIter     = ref.end();
        }
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
    IntVector next;
    // Advance the iterator. It is guaranteed not to be at the end yet.
    ++CurrentIter;

    // If we are out of particles in the current cell, search
    // through adjacent ones until we find another populated cell.
    while (CurrentIter == EndIter)
    {

        // Increment x
        ++CurrentCell[2];

        if (CurrentCell[2] > HighBounds[2]) // Carry over if x is out of bounds
        {
            ++CurrentCell[1];
            CurrentCell[2] = LowBounds[2];

            if (CurrentCell[1] > HighBounds[1]) // Carry over if y is out of bounds
            {
                ++CurrentCell[0];
                CurrentCell[1] = LowBounds[1];

                if (CurrentCell[0] > HighBounds[0]) // End if z is out of bounds.
                {
                    break;
                }
            }
        }

        next = this->Cell + this->CurrentCell;

        std::set<int>& ref = Parent.Cells[next[0]][next[1]][next[2]];

        if (!ref.empty())   // Only get the iterator if they are not empty.
        {
            CurrentIter = ref.begin();
            EndIter     = ref.end();
        }
    }

    return *this;
}

/**
 * @return The current value the iterator is focused on.
 */
int
CellList::NeighborCellIterator::operator*()
{
    return *CurrentIter;
}

/**
 * @return Returns true if the iterator can be safely dereferenced.
 */
bool
CellList::NeighborCellIterator::HasNext() const
{
    // The iterator is valid if the current adjacent cell is within bounds
    // and the iterator is not exhausted.
    return this->CurrentIter != this->EndIter;
}

void
CellList::Resize(const int ps, const int xs, const int ys, const int zs)
{
    // Resize to zero to clear out the arrays. Cheating? Perhaps, but it is
    // easy.
    Cells.resize(boost::extents[0][0][0]);
    Index.resize(boost::extents[0]);

    // Resize to the desired size.
    Cells.resize(boost::extents[xs][ys][zs]);
    Index.resize(boost::extents[ps]);

    // Nullify index so that we do not hide errors.
    IntVector* start            = Index.data();
    IntVector const * const end = start + Index.num_elements();

    for (; start != end; ++start)
    {
        (*start)[0] = (*start)[1] = (*start)[2] = -1;
    }

    Particles = ps;
    Size[0]   = xs;
    Size[1]   = ys;
    Size[2]   = zs;
}

void CellList::Clear()
{
    // Resize to zero to clear out the arrays. Cheating? Perhaps, but it is
    // easy.
    Cells.resize(boost::extents[0][0][0]);
    Index.resize(boost::extents[0]);

    // Resize to the desired size.
    Cells.resize(boost::extents[Size[0]][Size[1]][Size[2]]);
    Index.resize(boost::extents[Particles]);

    // Nullify index so that we do not hide errors.
    IntVector* start            = Index.data();
    IntVector const * const end = start + Index.num_elements();

    for (; start != end; ++start)
    {
        (*start)[0] = (*start)[1] = (*start)[2] = -1;
    }
}
