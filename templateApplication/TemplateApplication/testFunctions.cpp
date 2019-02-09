#include <cmath>

#include <TemplateApplication/testFunctions.h>

double performCalculation1(
    double input1, int input2 )
{
    return performCalculation2( input1, input2 ) - std::pow( input1,input2 );

}

double performCalculation2(
     double input1, double input2 )
{
    return input1 + input2;
}
