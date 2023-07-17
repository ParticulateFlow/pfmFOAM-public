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
    Foam::deviationPropagation::noDeviationPropagators

Description
    Empty class to store no deviation propagators.

\*---------------------------------------------------------------------------*/

#ifndef noDeviationPropagators_H
#define noDeviationPropagators_H

#include "deviationPropagators.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class singlePhaseDeviationPropagators Declaration
\*---------------------------------------------------------------------------*/

class noDeviationPropagators:
    public deviationPropagators
{
protected:


public:

    //- Runtime type information
    TypeName("noDeviationPropagators");


    // Constructors

        //- Construct from components
        noDeviationPropagators
        (
            const dictionary& dict,
            dataBase& base
        );

    // Destructor

        ~noDeviationPropagators();

    // Member Functions

    label readDeviationPropagators(fileNameList dataBases) {}

    tensorList &Kuu_boundary(label refState, label receiverID) {}

    tensorList &Kuu_internal(label refState, label receiverID) {}

    tensor Kuu_integrated(label refState, label cellI) {}

private:


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

