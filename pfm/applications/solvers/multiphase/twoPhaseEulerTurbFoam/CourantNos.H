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
    scalar CoNumP1(0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue());
    Info<< "Max Courant Number of phase 1 = " << CoNumP1 << endl;
    
    sumPhi = fvc::surfaceSum(mag(phi2))().primitiveField();
    scalar CoNumP2(0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue());
    Info<< "Max Courant Number of phase 2 = " << CoNumP2 << endl;
    

    CoNum = max(max(CoNum, CoNumP1),CoNumP2);
}
