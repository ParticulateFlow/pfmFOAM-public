#include "CourantNo.H"

/*
{
    scalar UrCoNum = 0.5*gMax
    (
        fvc::surfaceSum(mag(phi1 - phi2))().primitiveField()/mesh.V().field()
    )*runTime.deltaTValue();

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    CoNum = max(CoNum, UrCoNum);
}
*/
{
    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi1))().primitiveField()
    );
    
    sumPhi = max
    (
        sumPhi,
        fvc::surfaceSum(mag(phi2))().primitiveField()
    );

    scalar CoNumP(0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue());
    Info<< "Max Phase Courant Number = " << CoNumP << endl;

    CoNum = max(CoNum, CoNumP);
}
