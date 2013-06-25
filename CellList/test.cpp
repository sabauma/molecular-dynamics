#include <stdio.h>
#include <vector>

#include "CellListFast.h"

int main(void)
{
    CellList cells(100, 100, 100, 100);
    IntVector low(-1, -1, -1);
    IntVector high(2, 2, 2);

    cells.SetParticle(5, 5, 5, 5);

    for (int i = 0; i < 100; i += 2)
    {
        cells.SetParticle(i, i, i, i);
    }

    for (int i = 0; i < 50; ++i)
    {

        printf("Around particle %d\n", i);
        for (CellList::NeighborCellIterator iter =
                 cells.BeginNeighborCell(i*2, low, high);
             iter.HasNext(); ++iter)
        {
            printf("    Paricle %d\n", *iter);
        }
    }
}
