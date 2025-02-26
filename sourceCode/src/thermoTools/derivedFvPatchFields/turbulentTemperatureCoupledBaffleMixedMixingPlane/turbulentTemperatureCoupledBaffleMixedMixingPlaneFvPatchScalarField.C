/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()), //default method (fluidThermo)
    mappedPatchFieldBase<scalar>
    (
       mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this
    ),	
    TnbrName_("undefined-Tnbr"),
    fluidFilter_()
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        ptf
    ),
    TnbrName_(ptf.TnbrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    thicknessLayer_(ptf.thicknessLayer_.clone(p.patch())),
    kappaLayers_(ptf.kappaLayers_),
    kappaLayer_(ptf.kappaLayer_.clone(p.patch())),
    fluidFilter_(ptf.fluidFilter_)
{}

turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(p, iF),
        *this,
        dict
    ),
    TnbrName_(dict.get<word>("Tnbr")),
    fluidFilter_(dict.lookupOrDefault<bool>("fluidFilter","fluidFilter"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.readIfPresent("thicknessLayers", thicknessLayers_))
    {
        dict.readEntry("kappaLayers", kappaLayers_);
    }
    // Read single additional PatchFunction1
    thicknessLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "thicknessLayer",
        dict
    );
    kappaLayer_ = PatchFunction1<scalar>::NewIfPresent
    (
        p.patch(),
        "kappaLayer",
        dict
    );
     
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1.0;
    }
}


turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(wtcsf, iF),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), iF),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    thicknessLayer_(wtcsf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(wtcsf.kappaLayers_),
    kappaLayer_(wtcsf.kappaLayer_.clone(patch().patch())),
    fluidFilter_(wtcsf.fluidFilter_)
{}

turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::
turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
(
    const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField& wtcsf
)
:
    mixedFvPatchScalarField(wtcsf),
    temperatureCoupledBase(patch(), wtcsf),
    mappedPatchFieldBase<scalar>
    (
        mappedPatchFieldBase<scalar>::mapper(patch(), wtcsf.internalField()),
        *this,
        wtcsf
    ),
    TnbrName_(wtcsf.TnbrName_),
    thicknessLayers_(wtcsf.thicknessLayers_),
    thicknessLayer_(wtcsf.thicknessLayer_.clone(patch().patch())),
    kappaLayers_(wtcsf.kappaLayers_),
    kappaLayer_(wtcsf.kappaLayer_.clone(patch().patch())),
    fluidFilter_(wtcsf.fluidFilter_)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::autoMap(mapper);
    temperatureCoupledBase::autoMap(mapper);
    //mappedPatchFieldBase<scalar>::autoMap(mapper);
    if (thicknessLayer_)
    {
        thicknessLayer_().autoMap(mapper);
        kappaLayer_().autoMap(mapper);
    }
}


void turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::rmap
(
    const fvPatchField<scalar>& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField& tiptf =
        refCast
        <
            const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
        >(ptf);

    temperatureCoupledBase::rmap(tiptf, addr);
    //mappedPatchFieldBase<scalar>::rmap(ptf, addr);
    if (thicknessLayer_)
    {
        thicknessLayer_().rmap(tiptf.thicknessLayer_(), addr);
        kappaLayer_().rmap(tiptf.kappaLayer_(), addr);
    }
}


tmp<Foam::scalarField>
turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::kappa
(
    const scalarField& Tp
) const
{
    // Get kappa from relevant thermo
    tmp<scalarField> tk(temperatureCoupledBase::kappa(Tp));

    // Optionally modify with explicit resistance
    if (thicknessLayer_ || thicknessLayers_.size())
    {
        scalarField KDelta(tk*patch().deltaCoeffs());

        // Harmonic averaging of kappa*deltaCoeffs
        {
            KDelta = 1.0/KDelta;
            if (thicknessLayer_)
            {
                const scalar t = db().time().timeOutputValue();
                KDelta +=
                    thicknessLayer_().value(t)
                   /kappaLayer_().value(t);
            }
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    KDelta += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
                }
            }
            KDelta = 1.0/KDelta;
        }

        // Update kappa from KDelta
        tk = KDelta/patch().deltaCoeffs();
    }

    return tk;
}

void turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
    	mappedPatchFieldBase<scalar>::mapper
        (
            patch(),
            this->internalField()
        );
    const scalarField& Tp = *this;
    const scalarField kappaTp(kappa(Tp));
    const tmp<scalarField> myKDelta = kappaTp*patch().deltaCoeffs();

    scalarField nbrIntFld;
    scalarField nbrKDelta;
 
    // Initialise nbrField and nbrPatch
    const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const label nbrPatchID = mpp.samplePolyPatch().index(); 
    const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];

    const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField&
    nbrField =
    refCast
    <
        const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
    >
    (

        nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
    );

    if (mpp.sameWorld())
    {
        // Same world so lookup
//        const auto& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
//        const label nbrPatchID = mpp.samplePolyPatch().index(); 
//	const auto& nbrPatch = nbrMesh.boundary()[nbrPatchID];
//
//        const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField&
//        nbrField =
//        refCast
//        <
//        const turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
//        >
//        (
//            nbrPatch.lookupPatchField<volScalarField>(TnbrName_)
//        );

        // Swap to obtain full local values of neighbour K*delta
	nbrIntFld = nbrField.patchInternalField();
        nbrKDelta = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        // Different world so use my region,patch. Distribution below will
        // do the reordering.
	nbrIntFld = patchInternalField();
        nbrKDelta = myKDelta.ref();
    }
    distribute(this->internalField().name() + "_value", nbrIntFld);
    distribute(this->internalField().name() + "_weights", nbrKDelta);


    //Add time loop above this loop (if update,else reflect)
    //Update if statement (if nbr region is solid, else if nbr region is fluid)
    if (nbrField.fluidFilter_)
    {
        if(debug)
        {
            //Mixing Plane alert
            Info<< "Beginning Mixing Plane" << endl;

            //Debug: Info out
            Pout << "nbrPatch size: " << nbrPatch.size() << endl;
            Pout << "nbrPatch name: " << nbrPatch.name() << endl;
        }

        const Time& runTime = this->patch().boundaryMesh().mesh().time();
            //Create different lists from all processors
            //Number of Processors
            label n = Pstream::nProcs();
            label proci = Pstream::myProcNo();
            //Create lists
            DynamicList<double> procxList, procyList, proczList, procTList, procIDList;
            if(nbrPatch.Cf().size() != 0)
            {
                for(label ID = 0; ID < nbrPatch.Cf().size(); ID++)
                {
                    procxList.append(nbrPatch.Cf()[ID][0]);
                    procyList.append(nbrPatch.Cf()[ID][1]);
                    proczList.append(nbrPatch.Cf()[ID][2]);
                    procIDList.append(ID);
                }
            }
            //If size = 0 do nothing else do something
            if(nbrIntFld.size() != 0)
            {
                for(label ID = 0; ID < nbrIntFld.size(); ID++)
                {
                    procTList.append(nbrIntFld[ID]);
                }
            }
            //Create size list
            List<double> procSizeList(n, Zero);
            double procSize = patch().size();

            if(debug)
            {
                //Printing list elements for each processor
                Pout << "Processor " << proci << ": " << nl;
                Pout << "xList " << procxList << nl;
                Pout << "yList " << procyList << nl;
                Pout << "zList " << proczList << nl;
                Pout << "TList " << procTList << nl;
                Pout << "IDList " << procIDList << nl;
                Pout << "SizeList " << procSizeList << nl;
                Info << nl << endl;
            }

            //Create a placeholder list
            List<List<double>> sharexList(n), shareyList(n), sharezList(n), shareTList(n), shareIDList(n);
            List<double> shareSizeList(n);
            //Transfer prociList into the proci'th element of placeholder list shareList.
            //Removes content of prociList.
            sharexList[proci].transfer(procxList);
            shareyList[proci].transfer(procyList);
            sharezList[proci].transfer(proczList);
            shareTList[proci].transfer(procTList);
            shareIDList[proci].transfer(procIDList);
            shareSizeList[proci] = procSize;

            //Gather the different shareLists of slave processors in master processor
            Pstream::gatherList(sharexList);
            Pstream::gatherList(shareyList);
            Pstream::gatherList(sharezList);
            Pstream::gatherList(shareTList);
            Pstream::gatherList(shareIDList);
            Pstream::gatherList(shareSizeList);
            //Share masters shareList to the rest of the processors
            Pstream::scatterList(sharexList);
            Pstream::scatterList(shareyList);
            Pstream::scatterList(sharezList);
            Pstream::scatterList(shareTList);
            Pstream::scatterList(shareIDList);
            Pstream::scatterList(shareSizeList);
            //Create a single list with all elements in
            forAll(sharexList, listID)
            {
                procxList.append(sharexList[listID]);
                procyList.append(shareyList[listID]);
                proczList.append(sharezList[listID]);
                procTList.append(shareTList[listID]);
                procIDList.append(shareIDList[listID]);
            }
            forAll(shareSizeList, listID)
            {
                procSizeList[listID] = shareSizeList[listID];
            }

            if(debug)
            {
                //Debug: print final output list
                Pout<< "size of resulting procxList: " << procxList.size() << nl
                    << "procxList in proc " << proci << ": " << procxList << nl;

                Pout<< "size of resulting procyList: " << procyList.size() << nl
                    << "procyList in proc " << proci << ": " << procyList << nl;

                Pout<< "size of resulting proczList: " << proczList.size() << nl
                    << "proczList in proc " << proci << ": " << proczList << nl;

                Pout<< "size of resulting procTList: " << procTList.size() << nl
                    << "procTList in proc " << proci << ": " << procTList << nl;

                Pout<< "size of resulting procIDList: " << procIDList.size() << nl
                    << "procIDList in proc " << proci << ": " << procIDList << nl;

                Pout<< "size of resulting procSizeList: " << procSizeList.size() << nl
                    << "procSizeList in proc " << proci << ": " << procSizeList << nl;

                Pout<< "END PARALLEL DEBUGGING" << endl;

                //--------------------------------------------------------------------------//

                //Debug: Info out
                Info<< "nbrPatch size: " << nbrPatch.size() << endl;
                Info<< "nbrPatch name: " << nbrPatch.name() << endl;
            }

            //Initialise data Input and Output vectors
            std::vector<dataInput> dataIn;
            std::vector<dataOutput> dataOut;

            for(double faceID = 0; faceID < procxList.size(); faceID++)
            {
                dataIn.push_back(
                            dataInput(
                                procxList[faceID],
                                procyList[faceID],
                                proczList[faceID],
                                procTList[faceID],
                                procIDList[faceID],
                                faceID
                                ));

                if(debug)
                {
                    //Debug display dataIn population
                    Info<< dataIn[faceID].xIn << ": "
                        << dataIn[faceID].yIn << ": "
                        << dataIn[faceID].zIn << ": "
                        << dataIn[faceID].TIn << ": "
                        << dataIn[faceID].procIndexIn << ": "
                        << dataIn[faceID].indexIn
                        << endl;
                }
            }

            if(debug)
            {
                //DEBUG FLAG UPDATE
                Info<< "DataInput Allocation Successful" << endl;
            }
            //Run mixing plane function:
            bool mixingPlaneResult = mixingPlane(dataIn, &dataOut);

            if(debug)
            {
                //DEBUG FLAG UPDATE
                for(int j = 0; j < dataOut.size(); j++)
                {
                    Info<< dataOut[j].xOut << ": "
                        << dataOut[j].yOut << ": "
                        << dataOut[j].zOut << ": "
                        << dataOut[j].TOut << ": "
                        << dataOut[j].procOrigIndexOut << ": "
                        << dataOut[j].origIndexOut << ": "
                        << dataOut[j].rOut << ": "
                        << dataOut[j].thetaOut << ": "
                        << dataOut[j].TnewOut
                        << endl;
                }
                Info<< "mixingPlane Function Successful" << endl;
            }
            //Check for failed mixingPlane
            if (!mixingPlaneResult)
            {
                FatalErrorInFunction
                        << "mixingPlane Function failed"
                        << exit(FatalError);
            }

            //---------Repopulate T values - Tout is processor patch sized--------//
            //Create matrix of solid sizes
            List<double> sumSizeList(n);
            double sumSize = 0.0;
            for(label q=0; q < n; q++)
            {
                for(label w=0; w<q; w++)
                {
                    sumSize += procSizeList[w];
                }
                sumSizeList[q] = sumSize;
                sumSize = 0.0;
            }
            //Create field equal to solid patch size on *this* processor
            tmp<scalarField> TOut(new scalarField(patch().size(), Zero));
            int sumRef = sumSizeList[proci];
            for( label r = 0; r < patch().size(); r++)
            {
                TOut.ref()[r] = dataOut[sumRef+r].TOut;
            }

            if(debug)
            {
                //DEBUG FLAG UPDATE
                Info<< "TOut Allocation Successful" << endl;
            }
            //---------------------------------------------------------------------//

            //Clean
            dataIn.clear();
            dataOut.clear();

            //Set KDelta for this patch
            tmp<scalarField> myKDelta = kappa(*this)*patch().deltaCoeffs();

            //Update Temperature
            this->refValue() = TOut.ref();
            this->refGrad() = Zero;
            this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta());

            if(debug)
            {
                //Mixing Plane alert
                Info<< "Ending Mixing Plane" << endl;
            }
    }
    else
    {
        //Standard update method for coupled boundary, applied to the fluid domain
        //Standard Alert alert
	Info << "Standard Direct Coupling" << endl;

	// Both sides agree on
	// - temperature : (myKDelta*fld + nbrKDelta*nbrFld)/(myKDelta+nbrKDelta)
	// - gradient    : (temperature-fld)*delta
	// We've got a degree of freedom in how to implement this in a mixed bc.
	// (what gradient, what fixedValue and mixing coefficient)
	// Two reasonable choices:
	// 1. specify above temperature on one side (preferentially the high side)
	//    and above gradient on the other. So this will switch between pure
	//    fixedvalue and pure fixedgradient
	// 2. specify gradient and temperature such that the equations are the
	//    same on both sides. This leads to the choice of
	//    - refGradient = zero gradient
	//    - refValue = neighbour value
	//    - mixFraction = nbrKDelta / (nbrKDelta + myKDelta())

    	this->refValue() = nbrIntFld;
    	this->refGrad() = Zero;
    	this->valueFraction() = nbrKDelta/(nbrKDelta + myKDelta());
	

    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << mpp.sampleRegion() << ':'
            << mpp.samplePatch() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }

    // Restore tag
    UPstream::msgType(oldTag);
}



void turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& m,
    const label iMatrix,
    const direction cmpt
)
{
    FatalErrorInFunction
        << "This BC does not support energy coupling "
        << "Use compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << abort(FatalError);
}


tmp<Field<scalar>>
turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::coeffs
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& coeffs,
    const label mat
) const
{
    FatalErrorInFunction
        << "This BC does not support energy coupling "
        << "Use compressible::turbulentTemperatureRadCoupledMixed "
        << "which has more functionalities and it can handle "
        << "the assemble coupled option for energy. "
        << abort(FatalError);

    return nullptr;
}


void turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntry("Tnbr", TnbrName_);
    if (thicknessLayer_)
    {
        thicknessLayer_().writeData(os);
        kappaLayer_().writeData(os);
    }
    if (thicknessLayers_.size())
    {
        thicknessLayers_.writeEntry("thicknessLayers", os);
        kappaLayers_.writeEntry("kappaLayers", os);
    }
    temperatureCoupledBase::write(os);
    mappedPatchFieldBase<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //

