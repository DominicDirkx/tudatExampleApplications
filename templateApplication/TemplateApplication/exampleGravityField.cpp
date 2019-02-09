/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120522    A. Ronse          First creation of code.
 *
 *    References
 *      Williams, Dr. David R., "Moon Fact Sheet", NASA (National Space Science Data Center),
 *         http://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html, last accessed: 22 May 2012
 *
 *    Notes
 *
 */

#include <Eigen/Core>

#include <TemplateApplication/exampleGravityField.h>

// Constructor
GravityFieldModel::GravityFieldModel( double gravitationalParameter )
{
    gravitationalParameter_ = gravitationalParameter;
}

// Function to calculate the potential
double GravityFieldModel::calculatePotential( Eigen::Vector3d relativePosition )
{
    double distance = relativePosition.norm( );
    return -gravitationalParameter_ / distance;
}

// Function to calculate the potential
Eigen::Vector3d GravityFieldModel::calculatePotentialGradient( Eigen::Vector3d relativePosition )
{
    double distance = relativePosition.norm( );
    return relativePosition * gravitationalParameter_ / std::pow( distance, 3.0 );
}
