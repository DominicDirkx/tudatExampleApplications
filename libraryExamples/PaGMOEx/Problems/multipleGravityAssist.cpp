/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/propagationPatchedConicFullProblem.h"

#include "multipleGravityAssist.h"

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::transfer_trajectories;
using namespace tudat;
using namespace pagmo;

MultipleGravityAssist::MultipleGravityAssist(std::vector< std::vector< double > > &bounds,
                                             std::vector< int > flybySequence,
                                             const bool useTripTime ):
    problemBounds_( bounds ), useTripTime_( useTripTime )
{
    tudat::spice_interface::loadStandardSpiceKernels( );

    // Specify required parameters
    // Specify the number of legs and type of legs.
    numberOfLegs_ = flybySequence.size( );
    legTypeVector_.resize( numberOfLegs_ );
    legTypeVector_[ 0 ] = mga_Departure;
    legTypeVector_[ numberOfLegs_ - 1 ] = capture;

    for(int i = 1; i < numberOfLegs_ - 1; i++){
        legTypeVector_[ i ] = mga_Swingby;
    }

    // Create the ephemeris, gravitational parameter, and minimum pericentre vector.
    ephemerisVector_.resize( numberOfLegs_ );
    gravitationalParameterVector_.resize( numberOfLegs_ );
    minimumPericenterRadii_.resize( numberOfLegs_ );

    bodyNames.resize( numberOfLegs_ );
    for(int i = 0; i < numberOfLegs_; i++)
    {
        switch(flybySequence[ i ])
        {
        case( 1 ):
            bodyNames[ i ] = "Mercury";
            minimumPericenterRadii_[ i ] = 2639.7E3;
            break;
        case( 2 ):
            bodyNames[ i ] = "Venus";
            minimumPericenterRadii_[ i ] = 6251.8E3;
            break;
        case( 3 ):
            bodyNames[ i ] = "Earth";
            minimumPericenterRadii_[ i ] = 6578.1E3;
            break;
        case( 4 ):
            bodyNames[ i ] = "Mars";
            minimumPericenterRadii_[ i ] = 3596.2E3;
            break;
        case( 5 ):
            bodyNames[ i ] = "Jupiter";
            minimumPericenterRadii_[ i ] = 72000.0E3;
            break;
        case( 6 ):
            bodyNames[ i ] = "Saturn";
            minimumPericenterRadii_[ i ] = 61000.0E3;
            break;
        case( 7 ):
            bodyNames[ i ] = "Uranus";
            minimumPericenterRadii_[ i ] = 26000.0E3;
            break;
        case( 8 ):
            bodyNames[ i ] = "Neptune";
            minimumPericenterRadii_[ i ] = 25000.0E3;
            break;
        case( 9 ):
            bodyNames[ i ] = "Pluto";
            minimumPericenterRadii_[ i ] = 1395.0E3;
            break;
        default:
            std::cerr<<"Planet in flyby sequence is not defined.";
        }
    }

    bodyMap = propagators::setupBodyMapFromEphemeridesForPatchedConicsTrajectory(
                "Sun", "Vehicle", bodyNames );

    for(int i = 0; i < numberOfLegs_; i++)
    {
        ephemerisVector_[ i ] = bodyMap.at( bodyNames.at( i ) )->getEphemeris( );
        gravitationalParameterVector_[ i ] = bodyMap.at( bodyNames.at( i ) )->getGravityFieldModel( )->getGravitationalParameter( );
    }
    // Create departure and capture variables.
    semiMajorAxes_.resize( 2 );
    eccentricities_.resize( 2 );
    semiMajorAxes_ << std::numeric_limits< double >::infinity( ), 1.0895e8 / 0.02;
    eccentricities_ << 0., 0.98;

}
//! Descriptive name of the problem
std::string MultipleGravityAssist::get_name() const {
    return "MGA transfer trajectory";
}

//! Get bounds
std::pair<std::vector<double>, std::vector<double> > MultipleGravityAssist::get_bounds() const {

    return { problemBounds_[0], problemBounds_[1] };
}

//! Implementation of the fitness function (return delta-v)
std::vector<double> MultipleGravityAssist::fitness( const std::vector<double> &xv ) const{
    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;

    // Create variable vector.
    Eigen::VectorXd variableVector ( numberOfLegs_ + 1 );

    double TOF = 0;
    for(int i = 0; i < numberOfLegs_ ; i++){
        variableVector[ i ] = xv[ i ];
        if( i > 0 ){
            TOF += xv[i];
        }
    }
    variableVector[ numberOfLegs_ ] = 1;//dummy
    variableVector *= physical_constants::JULIAN_DAY;

     //Create the trajectory problem.
//        Trajectory mgaTraj( numberOfLegs_, legTypeVector_, ephemerisVector_,
//                              gravitationalParameterVector_, variableVector, sunGravitationalParameter,
//                              minimumPericenterRadii_, semiMajorAxes_, eccentricities_, false );
    Trajectory mgaTraj = createTransferTrajectoryObject(
                bodyMap, bodyNames, "Sun", legTypeVector_,
                utilities::convertEigenVectorToStlVector( variableVector ),
                utilities::convertEigenVectorToStlVector( minimumPericenterRadii_ ),
                false, TUDAT_NAN, TUDAT_NAN, true,
                1.0895e8 / 0.02, 0.98 );

    // Start the deltaV vector.
    double resultingDeltaV;
    mgaTraj.calculateTrajectory( resultingDeltaV );

    if (std::isnan(resultingDeltaV))
    {
        resultingDeltaV = 1.0E10;
    }

    if ( useTripTime_ ){
        return { resultingDeltaV, TOF };
    }
    else {
        return { resultingDeltaV };
    }

}


