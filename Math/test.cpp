
#include <vector>
#include <math.h>
#include "BSpline.h"

const double PI = 3.14159265358979323;

int main(int, char**)
{
    std::vector<double> x(100, 0.0);
    std::vector<double> y(100, 0.0);
        
    for (int i = 0; i < 100; ++i)
    {
        x[i] = 2 * PI / 100 * i;
        y[i] = sin(x[i]);
    }

    BSpline spline(12, 4, 12 + 2 - 4);
    spline.SetData(x, y);

    for (int i = 0; i < 200; ++i)
    {
        double xi = 2 * PI / 100 * i;
        printf("%e, %e, %e\n", xi, sin(xi), spline.Evaluate(xi));
    }

    return 0;
}
