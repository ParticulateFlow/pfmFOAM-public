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
    Foam::deviationPropagators::referenceStates::singlePhaseDeviationPropagators

Description
    Class to store deviation propagators for single-phase flow.

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseDeviationPropagators_H
#define singlePhaseDeviationPropagators_H

#include "deviationPropagators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class singlePhaseDeviationPropagators Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseDeviationPropagators:
    public deviationPropagators
{
protected:

    dictionary propsDict_;

    List<List<tensorList>> Kuu_boundary_;

    List<List<tensorList>> Kuu_internal_;

    List<tensorList> Kuu_integrated_;

public:

    //- Runtime type information
    TypeName("singlePhaseDeviationPropagators");


    // Constructors

        //- Construct from components
        singlePhaseDeviationPropagators
        (
            const dictionary& dict,
            dataBase& base
        );

    // Destructor

        ~singlePhaseDeviationPropagators();

    // Member Functions

    void readDeviationPropagators(fileNameList dataBases);

    tensorList &Kuu_boundary(label refState, label receiverID);

    tensorList &Kuu_internal(label refState, label receiverID);

    tensor Kuu_integrated(label refState, label cellI);

private:


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

