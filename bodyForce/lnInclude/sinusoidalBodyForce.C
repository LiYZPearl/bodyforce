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

#include "sinusoidalBodyForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(sinusoidalBodyForce, 0);
addToRunTimeSelectionTable(bodyForcing, sinusoidalBodyForce, bodyForcing);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

sinusoidalBodyForce::sinusoidalBodyForce
(
    const fvMesh & mesh,
    const dictionary & dict
)
:
    bodyForcing(mesh, dict),

    period_( readScalar( bodyForcingDict_.lookup("period") ) ),

    Um_( readScalar( bodyForcingDict_.lookup("Um") ) ),

    offset_( readScalar( bodyForcingDict_.lookup("constantSlope") ) ),

    direction_( vector( bodyForcingDict_.lookup("forceDirection" ) ))
{
    direction_ /= Foam::mag(direction_);
}


sinusoidalBodyForce::~sinusoidalBodyForce()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector sinusoidalBodyForce::correct( )
{
    scalar omega = 2 * M_PI / period_;

    dimensionedVector forcing
    (
        "forcing",
        dimensionSet(0,1,-2,0,0,0,0),
        (
            Um_ * omega * Foam::cos( omega * mesh_.time().value() )
          + offset_
        ) * direction_
    );

    return forcing;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
