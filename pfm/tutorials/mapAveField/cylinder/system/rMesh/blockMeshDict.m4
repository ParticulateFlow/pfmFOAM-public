// Parametrized test case for a rectangle


//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926536)
m4_define(sqrt2, 1.41)

//Geometry
// x-length before center
m4_define(x1, 0)
// x-length after center
m4_define(x2, 25)
// y-length above/below center
m4_define(y, 5)
// z-length
m4_define(z1, 0)
m4_define(z2, 0.5)


//Grid points (integers!):
m4_define(xNumberOfCells, 250)
m4_define(yNumberOfCells, 100)
m4_define(zNumberOfCells, 1)
m4_define(rGrading, 1.0)
//m4_define(rGrading, 0.5)


/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  6
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version         2.0;
    format          ascii;
    class           dictionary;
    object          blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    (x1 -y z1)                          vlabel(V0)
    (x2 -y z1)                          vlabel(V1)
    (x2 y z1)                           vlabel(V2)
    (x1 y z1)                           vlabel(V3)
    (x1 -y z2)                          vlabel(V4)
    (x2 -y z2)                          vlabel(V5)
    (x2 y z2)                           vlabel(V6)
    (x1 y z2)                           vlabel(V7)
);

// Defining blocks:
blocks
(
    hex (0 1 2 3 4 5 6 7) (xNumberOfCells yNumberOfCells zNumberOfCells) simpleGrading (1 1 1)
);

edges
(
);

// Defining patches:
boundary
(
    leftPatch
    {
        type patch;
        faces
        (
            (V0 V4 V7 V3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (V1 V5 V6 V2)
        );
    }

    topbottom
    {
        type wall;
        faces
        (
            (V0 V1 V5 V4)
            (V3 V2 V6 V7)
        );
    }


    front
    {
        type empty;
        //neighbourPatch  back;
        faces
        (
            (V0 V1 V2 V3)
        );
    }


    back
    {
        type empty;
        //neighbourPatch  front;
        faces
        (
            (V4 V5 V6 V7)
        );
    }

);

mergePatchPairs
(
);

// ************************************************************************* //
