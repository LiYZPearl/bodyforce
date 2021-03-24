/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "stokes2ndBodyForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(stokes2ndBodyForce, 0);
addToRunTimeSelectionTable(bodyForcing, stokes2ndBodyForce, bodyForcing);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stokes2ndBodyForce::stokes2ndBodyForce
(
    const fvMesh & mesh,
    const dictionary & dict
)
:
    bodyForcing(mesh, dict),

    period_( readScalar( bodyForcingDict_.lookup("period") ) ),

    Um1_( readScalar( bodyForcingDict_.lookup("Um1") ) ),

    Um2_( readScalar( bodyForcingDict_.lookup("Um2") ) ),

    offset_( readScalar( bodyForcingDict_.lookup("constantSlope") ) ),

    direction_( vector( bodyForcingDict_.lookup("forceDirection" ) ))
{
    direction_ /= Foam::mag(direction_);

    Info << period_ << tab << Um1_ << tab << offset_ << tab << direction_ << endl;
}


stokes2ndBodyForce::~stokes2ndBodyForce()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector stokes2ndBodyForce::correct( )
{
    scalar omega = 2 * M_PI / period_;

    dimensionedVector forcing
    (
        "forcing",
        dimensionSet(0,1,-2,0,0,0,0),
        (
                  Um1_ * omega * Foam::cos( omega * mesh_.time().value() )
          + 2.0 * Um2_ * omega * Foam::sin( omega * mesh_.time().value() )
          + offset_
        ) * direction_
    );

    return forcing;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
