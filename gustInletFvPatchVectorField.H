/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::gustInletFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Implementation started from rampedParabolicInletFvPatchVectorField.

    Provides a parabolic inlet profile (according to the parameters specified
    by Turek and Hron 2006. Ramping over time is applied with a 
    (1-cos(pi/2*t))/2 function. 

SeeAlso
    Foam::fixedValueFvPatchVectorField

SourceFiles
    gustInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef gustInletFvPatchVectorField_H
#define gustInletFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class movingWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class gustInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data
        vector      refValue_;

        //- gust description
        vector      gustDirection_;
        //scalarList  gustAmplitudes_; // "Not implemented" error on read?!?
        //scalarList  gustFrequencies_; // "Not implemented" error on read?!?
        scalar  gustAmplitude_;
        scalar  gustFrequency_;
            
        //- line source location/extent
        vector      r0_;
        vector      avec_;
        scalar      sourceRadius_;


    // Private Member Functions

        //- Return current scale
        scalar currentScale() const;

public:

    //- Runtime type information
    TypeName("gustInlet");


    // Constructors

        //- Construct from patch and internal field
        gustInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        gustInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given gustInletFvPatchVectorField
        //  onto a new patch
        gustInletFvPatchVectorField
        (
            const gustInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        gustInletFvPatchVectorField
        (
            const gustInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new gustInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        gustInletFvPatchVectorField
        (
            const gustInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new gustInletFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
