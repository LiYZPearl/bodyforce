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

#include "bodyForcing.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(bodyForcing, 0);
defineRunTimeSelectionTable(bodyForcing, bodyForcing);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bodyForcing::bodyForcing
(
    const fvMesh & mesh,
    const dictionary & dict
)
:
    mesh_(mesh),

    bodyForcingDict_
    (
        dict
    )
{

}


bodyForcing::~bodyForcing()
{}


autoPtr<bodyForcing> bodyForcing::New
(
    const fvMesh & mesh,
    const dictionary & dict
)
{
    word formulation;

    dict.lookup("forceModel") >> formulation;

    bodyForcingConstructorTable::iterator cstrIter =
            bodyForcingConstructorTablePtr_->find(formulation);

    if (cstrIter == bodyForcingConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "bodyForcing::New(const fvMesh &, const word &)"
        )   << "Unknown body force formulation: " << formulation
            << endl << endl
            << "Valid methods are :" << endl
            << bodyForcingConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<bodyForcing>(cstrIter()( mesh, dict ));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
