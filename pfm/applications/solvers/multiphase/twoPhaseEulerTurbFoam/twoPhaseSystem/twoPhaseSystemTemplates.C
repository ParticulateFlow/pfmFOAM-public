/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "BlendedInterfacialModel.H"
#include "dragModel.H"
#include "driftVelocityModel.H"
#include "interPhaseForceModel.H"
#include "driftTemperatureModel.H"
#include "heatTransferModel.H"
#include "virtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class modelType>
const modelType& twoPhaseSystem::lookupSubModel
(
    const phasePair& key
) const
{
    return
        mesh().lookupObject<modelType>
        (
            IOobject::groupName(modelType::typeName, key.name())
        );
}


template<>
inline const dragModel& twoPhaseSystem::lookupSubModel<dragModel>
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return drag_->phaseModel(dispersed);
}

template<>
inline const heatTransferModel& twoPhaseSystem::lookupSubModel<heatTransferModel>
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return heatTransfer_->phaseModel(dispersed);
}

template<>
inline const driftVelocityModel& twoPhaseSystem::lookupSubModel<driftVelocityModel>
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return driftVelocity_->phaseModel(dispersed);
}
    
template<>
inline const interPhaseForceModel& twoPhaseSystem::lookupSubModel<interPhaseForceModel>
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return interPhaseForce_->phaseModel(dispersed);
}

template<>
inline const driftTemperatureModel& twoPhaseSystem::lookupSubModel<driftTemperatureModel>
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return driftTemperature_->phaseModel(dispersed);
}
    
template<>
inline const virtualMassModel& twoPhaseSystem::lookupSubModel<virtualMassModel>
(
    const phaseModel& dispersed,
    const phaseModel& continuous
) const
{
    return virtualMass_->phaseModel(dispersed);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
