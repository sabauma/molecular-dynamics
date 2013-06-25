
#include "rapidjson/prettywriter.h"
#include "rapidjson/filestream.h"
#include <stdio.h>
#include "InfoStruct.h"

int main(void)
{
    using namespace rapidjson;
    InfoStruct is = InfoStruct::ParseInfoStruct("./Example.json");

    printf("%0.5g\n", is.ConeAngle);
    printf("%0.5g\n", is.Temperature);
    printf("%0.5g\n", is.WallTemp);
    printf("%0.5g\n", is.R0);
    printf("%0.5g\n", is.Pd);
    printf("%0.5g\n", is.Frequency);
    printf("%0.5g\n", is.Tini);
    printf("%0.5g\n", is.Rini);
    printf("%0.5g\n", is.Vini);
    printf("%0.5g\n", is.Radi);
    printf("%0.5g\n", is.Riso);

    for (unsigned int i = 0; i < is.AtomicParts.size(); ++i)
    {
        printf("***");
        printf("%0.3g\n", is.AtomicParts[i]);
    }

    return 0;
}
