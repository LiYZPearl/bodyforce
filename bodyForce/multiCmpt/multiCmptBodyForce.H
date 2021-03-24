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
    Foam::multiCmptBodyForce

Description
    Returns a multi component motion. Includes the option of an additional
    constant body force.

SourceFiles
    multiCmptBodyForce.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef multiCmptBodyForce_H
#define multiCmptBodyForce_H

#include "bodyForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class multiCmptBodyForce
:
    public bodyForcing
{

protected:

    // Protected data
        scalarField period_;
        vectorField Um_;
        vector constant_;

    // Protected member functions

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiCmptBodyForce(const multiCmptBodyForce&);

        //- Disallow default bitwise assignment
        void operator=(const multiCmptBodyForce&);

public:

    //- Runtime type information
    TypeName("multiCmptBodyForce");

    // Constructors

        //- Construct from components
        multiCmptBodyForce
        (
            const fvMesh &,
            const dictionary &
        );

    // Destructor

        virtual ~multiCmptBodyForce();

    // Member Functions

        //- Compute bed load transport vector
        virtual dimensionedVector correct( );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
