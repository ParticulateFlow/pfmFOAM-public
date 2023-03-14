/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Contributing authors:
    Thomas Lichtenegger
    Copyright (C) 2015- Johannes Kepler University, Linz
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
    Foam::deviationResponse::referenceStates::generalReferenceStates

Description
    Class to store reference states for arbitrary flows.

\*---------------------------------------------------------------------------*/

#ifndef generalReferenceStates_H
#define generalReferenceStates_H

#include "referenceStates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class generalReferenceStates Declaration
\*---------------------------------------------------------------------------*/

class generalReferenceStates:
    public referenceStates
{
protected:

    bool verbose_;

public:

    //- Runtime type information
    TypeName("generalReferenceStates");


    // Constructors

        //- Construct from components
        generalReferenceStates
        (
            const dictionary& dict,
            dataBase& base
        );


    // Destructor

        ~generalReferenceStates();


    // Member Functions

    label readReferenceStates(wordList dataBases);

    const volScalarField& exportVolScalarField(label, label);

    const volVectorField& exportVolVectorField(label, label);

    const volScalarField& exportVolScalarEvolvedField(label, label);

    const volVectorField& exportVolVectorEvolvedField(label, label);

private:

    label maxNumRefStates_;

    List<PtrList<volScalarField>> volScalarRefStateList_;

    List<PtrList<volVectorField>> volVectorRefStateList_;

    List<PtrList<volScalarField>> volScalarEvolvedRefStateList_;

    List<PtrList<volVectorField>> volVectorEvolvedRefStateList_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

