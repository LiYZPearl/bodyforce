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

#include "multiCmptBodyForce.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(multiCmptBodyForce, 0);
addToRunTimeSelectionTable(bodyForcing, multiCmptBodyForce, bodyForcing);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

multiCmptBodyForce::multiCmptBodyForce
(
    const fvMesh & mesh,
    const dictionary & dict
)
:
    bodyForcing(mesh, dict),

    period_( 0, 0.0 ),

    Um_( 0, vector::zero ),

    constant_( bodyForcingDict_.lookup("constantSlope") )
{
    scalarField pp( bodyForcingDict_.lookup("period") );
    vectorField uu( bodyForcingDict_.lookup("Um") );

    period_.setSize(pp.size());
    period_ = pp;

    Um_.setSize(uu.size());
    Um_ = uu;
}


multiCmptBodyForce::~multiCmptBodyForce()
{}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector multiCmptBodyForce::correct( )
{
    dimensionedVector forcing
    (
        "forcing",
        dimensionSet(0,1,-2,0,0,0,0),
        constant_
    );

    forAll(period_, cmpt)
    {
        scalar omega = 2 * M_PI / period_[cmpt];
        forcing.value() += Um_[cmpt] * omega * Foam::cos( omega * mesh_.time().value() );
    }

    return forcing;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
