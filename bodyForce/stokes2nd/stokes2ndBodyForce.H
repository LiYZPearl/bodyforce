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
    Foam::stokes2ndBodyForce

Description
    Adds a simple second order correction to the sinusoidal body force.
    Also includes the option for a constant offset, i.e. current.

SourceFiles
    stokes2ndBodyForce.C

Author
    Niels Gjøl Jacobsen, Technical University of Denmark.
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef stokes2ndBodyForce_H
#define stokes2ndBodyForce_H

#include "bodyForcing.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveTheory Declaration
\*---------------------------------------------------------------------------*/

class stokes2ndBodyForce
:
    public bodyForcing
{

protected:

    // Protected data
        scalar period_;
        scalar Um1_;
        scalar Um2_;
        scalar offset_;

        vector direction_;

    // Protected member functions

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        stokes2ndBodyForce(const stokes2ndBodyForce&);

        //- Disallow default bitwise assignment
        void operator=(const stokes2ndBodyForce&);

public:

    //- Runtime type information
    TypeName("stokes2ndBodyForce");

    // Constructors

        //- Construct from components
        stokes2ndBodyForce
        (
            const fvMesh &,
            const dictionary &
        );

    // Destructor

        virtual ~stokes2ndBodyForce();

    // Member Functions

        //- Compute bed load transport vector
        virtual dimensionedVector correct( );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
