/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <limits>

#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h>
#include <Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h>

#include "Pagmo2Example/Problems/earthMarsTransfer.h"

using namespace pagmo;

//EarthMarsTransfer::EarthMarsTransfer(
//    const std::vector< std::vector< double > > problemBounds ) :
//    problemBounds_( problemBounds )
//{ }

//! Descriptive name of the problem
std::string EarthMarsTransfer::get_name() const {
    return "Multi-revolution Lambert Earth-Mars transfer trajectory";
}

//! Get bounds
std::pair<vector_double, vector_double> EarthMarsTransfer::get_bounds() const {
   std::vector< std::vector< double > > bounds( 2, std::vector< double >( 2, 0.0 ) );

    // Search between 2020 and 2025 for flight duration between 200
    // and 1000 days.
    bounds[ 0 ][ 0 ] = 2458849.5;
    bounds[ 1 ][ 0 ] = 2460676.5;
    bounds[ 0 ][ 1 ] = 200;
    bounds[ 1 ][ 1 ] = 1000;

//    return {problemBounds_[0], problemBounds_[1]};
    return {bounds[0], bounds[1]};
}

//! Implementation of the fitness function (return delta-v)
vector_double EarthMarsTransfer::fitness( const vector_double &xv ) const{

    using tudat::mission_segments::MultiRevolutionLambertTargeterIzzo;

    vector_double f;

    f.resize( 1 );

    // Gravitational parameter of the Sun
    double mu = 1.32712440018e+20;

    // Set initial and final position as those of Earth and Mars at
    // departure and arrival respectively.
    StateType initialState = getPlanetPosition( xv[0], "Earth");
    StateType finalState   = getPlanetPosition( xv[0] + xv[1], "Mars" );

    MultiRevolutionLambertTargeterIzzo lambertTargeter( initialState.segment(0,3),
        finalState.segment(0,3), xv[1]*86400, mu );
    double deltaV = std::numeric_limits<double>::infinity();
    unsigned int maxrev = lambertTargeter.getMaximumNumberOfRevolutions( );

    // Go through all multi-revolution solutions and select the one
    // with the lowest delta-V
    for( unsigned int i = 0; i <= maxrev; ++i){
	lambertTargeter.computeForRevolutionsAndBranch( i, false );
	deltaV = std::min( deltaV, ( initialState.segment(3,3)
		- lambertTargeter.getInertialVelocityAtDeparture( )).norm() +
	    + ( finalState.segment(3,3)
		- lambertTargeter.getInertialVelocityAtArrival( )).norm());
    }
    f[0] = deltaV;

    return f;
}

//! Function to obtain position of Earth and Mars
StateType EarthMarsTransfer::getPlanetPosition( const double date,
                                                const std::string planetName ) const {
    using tudat::orbital_element_conversions::convertKeplerianToCartesianElements;
    using tudat::orbital_element_conversions::convertMeanAnomalyToEccentricAnomaly;
    using tudat::orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly;


    // Gravitational parameter of the Sun
    double mu = 1.32712440018e+20;

    StateType stateKepl, stateCart;
    double n, jd0;
    if( planetName == "Earth" ){
        n   = 1.991e-07;
        jd0 = 2454000.5;
        stateKepl << 1.4960e+11, 1.6717e-02, 0.0, 5.0198e+00, 3.0614e+00, 4.4961e+00;
    } else {
        n   = 1.0586e-07;
        jd0 = 2451545.0;
        stateKepl << 2.2794e+11, 9.3412e-02, 0.0, 5.8650e+00, 8.6531e-01, 5.7567e+00;
    }
    stateKepl( 5 ) = convertMeanAnomalyToEccentricAnomaly( stateKepl( 1 ),
        fmod( stateKepl( 5 ) + ( date - jd0 ) * 86400. * n, 2.*M_PI ) );
    stateKepl( 5 ) = convertEccentricAnomalyToTrueAnomaly( stateKepl( 5 ), stateKepl( 1 ) );
    stateCart = convertKeplerianToCartesianElements( stateKepl , mu );
    return stateCart;
}
