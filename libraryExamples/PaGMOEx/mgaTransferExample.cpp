/*    Copyright (c) 2010-2017, Delft University of Technology
*    All rigths reserved
*
*    This file is part of the Tudat. Redistribution and use in source and
*    binary forms, with or without modification, are permitted exclusively
*    under the terms of the Modified BSD license. You should have received
*    a copy of the license with this file. If not, please or visit:
*    http://tudat.tudelft.nl/LICENSE.
*/

#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>

#include "pagmo/algorithms/de1220.hpp"
#include "pagmo/algorithms/sade.hpp"
#include "pagmo/algorithms/simulated_annealing.hpp"
#include "pagmo/algorithms/cmaes.hpp"
#include "pagmo/algorithms/de.hpp"
#include "pagmo/algorithms/nsga2.hpp"

#include "Problems/multipleGravityAssist.h"
#include "Problems/applicationOutput.h"
#include "Problems/getAlgorithm.h"
#include "Problems/saveOptimizationResults.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"

//! Execute  main
int main( )
{
    //Set seed for reproducible results
    pagmo::random_device::set_seed( 123456789 );

    std::vector< std::pair< int, int > > bodyIndices =
    { { 3, 3 }, { 2, 3 }, { 3, 2 }, { 2, 4 }, { 3, 4 }, { 4, 4 }, { 4, 2 } };
    std::vector< std::string > bodyLetters =
    { "EE", "VE", "EV", "VM", "EM", "MM", "MV" };

    for( int j = 0; j < bodyLetters.size( ); j++ )
    {
        // We have five decision variables each with a lower and upper
        // bound, create a vector of vectors that will contain these.
        int numberOfParameters = 5;
        std::vector< std::vector< double > > bounds( 2, std::vector< double >( numberOfParameters, 0.0 ) );

        // Define search bounds: first parameter is start date, following parameters are leg durations
        bounds[ 0 ][ 0 ] = -15 * 365; //MJD2000
        bounds[ 1 ][ 0 ] = 15 * 365; //MJD2000
        bounds[ 0 ][ 1 ] = 50;
        bounds[ 1 ][ 1 ] = 500;
        bounds[ 0 ][ 2 ] = 50;

        bounds[ 1 ][ 2 ] = 500;
        bounds[ 0 ][ 3 ] = 50;

        bounds[ 1 ][ 3 ] = 500;
        bounds[ 0 ][ 4 ] = 100;
        bounds[ 1 ][ 4 ] = 1600;

        if( j == 1 || j == 3 )
        {
            bounds[ 0 ][ 2 ] = 150;
        }

        if( j == 0 )
        {
            bounds[ 0 ][ 3 ] = 150;
        }

        if( j == 5 )
        {
            bounds[ 0 ][ 3 ] = 200;
            bounds[ 1 ][ 3 ] = 800;
        }

        // Define the problem: EVEEJ flyby sequence
        std::vector< int > flybySequence;
        flybySequence.push_back( 3 );
        flybySequence.push_back( 2 );
        flybySequence.push_back( bodyIndices.at( j ).first );
        flybySequence.push_back( bodyIndices.at( j ).second );
        flybySequence.push_back( 5 );


        // Create object to compute the problem fitness
        problem prob{ MultipleGravityAssist( bounds, flybySequence, true ) };

        // Select NSGA2 algorithm for priblem
        algorithm algo{nsga2( )};

        // Create an island with 1000 individuals
        island isl{algo, prob, 100 };

        // Evolve for 512 generations
        for( int i = 0 ; i < 1000; i++ )
        {

            isl.evolve( );
            while( isl.status( ) != pagmo::evolve_status::idle &&
                   isl.status( ) != pagmo::evolve_status::idle_error )
            {
                isl.wait( );
            }

            isl.wait_check( ); // Raises errors

            // Write current iteration results to file
            printPopulationToFile( isl.get_population( ).get_x( ), "mo_mga_EV" + bodyLetters.at( j ) + "J_" +
                                   std::to_string( i ), false );
            printPopulationToFile( isl.get_population( ).get_f( ), "mo_mga_EV" + bodyLetters.at( j ) + "J_" +
                                   std::to_string( i ), true );
            std::cout<<i<<std::endl;
        }
        double filterValue = 20.0E3;
        std::vector< std::vector< double > > finalParameters = isl.get_population( ).get_x( );
        std::vector< std::vector< double > > finalFitness = isl.get_population( ).get_f( );

        std::vector< std::vector< double > > filteredParameters;
        std::vector< std::vector< double > > filteredFitness;

        for( int i = 0; i < finalFitness.size( ); i++ )
        {
            if( finalFitness.at( i ).at( 0 ) < filterValue )
            {
                filteredParameters.push_back( finalParameters.at( i ) );
                filteredFitness.push_back( finalFitness.at( i ) );
            }
        }

        printPopulationToFile( filteredParameters, "mo_mga_EV" + bodyLetters.at( j ) + "J_filtered", false );
        printPopulationToFile( filteredFitness, "mo_mga_EV" + bodyLetters.at( j ) + "J_filtered", true );

    }
    return 0;

}
