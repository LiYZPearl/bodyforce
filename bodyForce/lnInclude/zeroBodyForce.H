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
    Foam::zeroBodyForce

Description
    A dummy implementation of the body force, where a 0 (zero) force is applied.

SourceFiles
    zeroBodyForce.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef zeroBodyForce_H
#define zeroBodyForce_H

#include "bodyForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class zeroBodyForce
:
    public bodyForcing
{

protected:

    // Protected data

    // Protected member functions

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        zeroBodyForce(const zeroBodyForce&);

        //- Disallow default bitwise assignment
        void operator=(const zeroBodyForce&);

public:

    //- Runtime type information
    TypeName("zeroBodyForce");

    // Constructors

        //- Construct from components
        zeroBodyForce
        (
            const fvMesh &,
            const dictionary &
        );

    // Destructor

        virtual ~zeroBodyForce();

    // Member Functions

        //- Compute bed load transport vector
        virtual dimensionedVector correct( );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
