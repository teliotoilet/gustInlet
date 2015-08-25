/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "gustInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
//#include "fvcMeshPhi.H"
#include "mathematicalConstants.H"
#include "steadyStateDdtScheme.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gustInletFvPatchVectorField::
gustInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(vector(1,0,0)),
    gustDirection_(vector(0,0,1)),
    //gustAmplitudes_(0),
    //gustFrequencies_(0),
    gustAmplitude_(0.0),
    gustFrequency_(0.0),
    r0_(vector::zero),
    avec_(vector(0,1,0)),
    sourceRadius_(1.0)
{}


Foam::gustInletFvPatchVectorField::
gustInletFvPatchVectorField
(
    const gustInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    refValue_(ptf.refValue_),
    gustDirection_(ptf.gustDirection_),
    //gustAmplitudes_(ptf.gustAmplitudes_),
    //gustFrequencies_(ptf.gustFrequencies_),
    gustAmplitude_(ptf.gustAmplitude_),
    gustFrequency_(ptf.gustFrequency_),
    r0_(ptf.r0_),
    avec_(ptf.avec_),
    sourceRadius_(ptf.sourceRadius_)
{}


Foam::gustInletFvPatchVectorField::
gustInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(dict.lookupOrDefault<vector>("Uinf",vector(1,0,0))),
    gustDirection_(dict.lookupOrDefault<vector>("direction",vector(0,0,1))),
//    gustAmplitudes_(dict.lookup("amplitude"))//,
//    gustFrequencies_(dict.lookup("frequency")), // Hz
    gustAmplitude_(dict.lookupOrDefault("amplitude",0.0)),
    gustFrequency_(dict.lookupOrDefault("frequency",0.0)), // Hz
    r0_(dict.lookupOrDefault<vector>("sourceLocation",vector(0,0,0))),
    avec_(dict.lookupOrDefault<vector>("sourceDirection",vector(0,1,0))),
    sourceRadius_(dict.lookupOrDefault<scalar>("sourceRadius",1.0))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::gustInletFvPatchVectorField::
gustInletFvPatchVectorField
(
    const gustInletFvPatchVectorField& pvf
)
:
    fixedValueFvPatchVectorField(pvf),
    refValue_(pvf.refValue_),
    gustDirection_(pvf.gustDirection_),
    //gustAmplitudes_(pvf.gustAmplitudes_),
    //gustFrequencies_(pvf.gustFrequencies_),
    gustAmplitude_(pvf.gustAmplitude_),
    gustFrequency_(pvf.gustFrequency_),
    r0_(pvf.r0_),
    avec_(pvf.avec_),
    sourceRadius_(pvf.sourceRadius_)
{}


Foam::gustInletFvPatchVectorField::
gustInletFvPatchVectorField
(
    const gustInletFvPatchVectorField& pvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pvf, iF),
    refValue_(pvf.refValue_),
    gustDirection_(pvf.gustDirection_),
    //gustAmplitudes_(pvf.gustAmplitudes_),
    //gustFrequencies_(pvf.gustFrequencies_),
    gustAmplitude_(pvf.gustAmplitude_),
    gustFrequency_(pvf.gustFrequency_),
    r0_(pvf.r0_),
    avec_(pvf.avec_),
    sourceRadius_(pvf.sourceRadius_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gustInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField::operator=( refValue_ );

    word ddtScheme
    (
        this->dimensionedInternalField().mesh()
        .ddtScheme(this->dimensionedInternalField().name())
    );

    // only add gust for unsteady case!
    if (ddtScheme != fv::steadyStateDdtScheme<scalar>::typeName)
    {
        const fvPatch& p = patch();
        const polyPatch& pp = p.patch();
        //const vectorField xc = pp.faceCentres();

        const Foam::scalar t = this->db().time().timeOutputValue();

        const vectorField rvec( pp.faceCentres() - r0_ );
        const vectorField bvec( rvec - (rvec & avec_)*avec_ );
        const scalarField expfn( Foam::exp(-(bvec & bvec)/(2*sourceRadius_*sourceRadius_)) );

        scalar scaling(0.0);
        scalar amplitude(0.0);
        scalar radPerSec(0.0);

    //    forAll(gustFrequencies_, mode)
    //    {
            amplitude = gustAmplitude_; //gustAmplitudes_[mode];
            radPerSec = 2*constant::mathematical::pi * gustFrequency_; //gustFrequencies_[mode];
            scaling = amplitude * sin( radPerSec*t );
            Info<< "gustInlet "// mode " << mode
                << " : amplitude=" << amplitude << " m/s"
                << ", frequency=" << radPerSec << " rad/s"
                << ", current scaling=" << scaling << endl;

            vectorField::operator+=( scaling * expfn * gustDirection_ );
    //    }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::gustInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Uinf") << refValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction") << gustDirection_ << token::END_STATEMENT << nl;
    //os.writeKeyword("amplitude") << gustAmplitudes_ << token::END_STATEMENT << nl;
    //os.writeKeyword("frequency") << gustFrequencies_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude") << gustAmplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << gustFrequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("sourceLocation") << r0_ << token::END_STATEMENT << nl;
    os.writeKeyword("sourceDirection") << avec_ << token::END_STATEMENT << nl;
    os.writeKeyword("sourceRadius") << sourceRadius_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        gustInletFvPatchVectorField
    );
}

// ************************************************************************* //
