/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::compressible::
        turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    on back-to-back baffles. Optional thin thermal layer resistances can be
    specified through thicknessLayers and kappaLayers entries.

    Specifies gradient and temperature such that the equations are the same
    on both sides:
      - refGradient = zero gradient
      - refValue = neighbour value
      - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    where KDelta is heat-transfer coefficient K * deltaCoeffs

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class temperatureCoupledBase.

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivities per layer [W/m/K] | no |
        kappaMethod  | inherited from temperatureCoupledBase | inherited |
        kappa        | inherited from temperatureCoupledBase | inherited |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureCoupledBaffleMixedMixingPlane;
        Tnbr            T;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        kappaMethod     lookup;
        kappa           kappa;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam::temperatureCoupledBase

SourceFiles
    turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField_H
#define turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "mappedPatchFieldBase.H"
#include "scalarField.H"
#include "turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "IOstream.H"
#include "cylindricalCS.H"
#include "mathematicalConstants.H"
#include "bndStructs.H"
#include "mixingPlane.H"
#include <vector>
#include <cmath>
#include <iterator>
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase,
    public mappedPatchFieldBase<scalar>
{
    // Private data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;
	autoPtr<PatchFunction1<scalar>> thicknessLayer_;

        //- Conductivity of layers
        scalarList kappaLayers_;
	autoPtr<PatchFunction1<scalar>> kappaLayer_;
    
        //- Fluid or Solid (for filtered boundary)
        bool fluidFilter_;

    // Private member functions

        //- Calculate coefficients for assembly matrix
        tmp<Field<scalar>> coeffs
        (
            fvMatrix<scalar>& matrix,
            const Field<scalar>&,
            const label
        ) const;

public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureCoupledBaffleMixedMixingPlane");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
        }

    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        //- Given patch temperature calculate corresponding K field. Override
        //- temperatureCoupledBase::kappa to includes effect of any
        //- explicit kappaThickness
        virtual tmp<scalarField> kappa(const scalarField& Tp) const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Manipulate matrix
        virtual void manipulateMatrix
        (
            fvMatrix<scalar>& m,
            const label iMatrix,
            const direction cmpt
        );


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
