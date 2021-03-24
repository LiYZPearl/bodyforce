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

Class
    Foam::bodyForce

Description
    An abstract base class for the computation of a body force to drive a
    flow in a periodic domain.

SourceFiles
    bodyForce.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef bodyForcing_H
#define bodyForcing_H

#include "fvMesh.H"
#include "volFields.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class bodyForcing
{

protected:

    // Protected data
        const fvMesh & mesh_;

        const dictionary & bodyForcingDict_;

    // Protected member functions

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        bodyForcing(const bodyForcing&);

        //- Disallow default bitwise assignment
        void operator=(const bodyForcing&);

public:

    //- Runtime type information
    TypeName("bodyForcing");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            bodyForcing,
            bodyForcing,
            (
                const fvMesh & mesh,
                const dictionary & dict
            ),
            ( mesh, dict )
        );

    // Constructors

        //- Construct from components
        bodyForcing
        (
            const fvMesh &,
            const dictionary &
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<bodyForcing> New
        (
            const fvMesh &,
            const dictionary &
        );

    // Destructor

        virtual ~bodyForcing();

    // Member Functions

        //- Return the sand slide modified bed level change
        virtual dimensionedVector correct(  ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
