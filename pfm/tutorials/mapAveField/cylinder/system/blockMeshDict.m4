// Parametrized test case for a cylinder geometry


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
// diameter of cylinder
m4_define(D, 1)
// outer diameter
m4_define(D2, 2)
// x-length before center
m4_define(x1, -5)
// x-length after center
m4_define(x2, 25)
// y-length above/below center
m4_define(y, 5)
// z-length
m4_define(z1, 0)
m4_define(z2, 0.5)








//Grid points (integers!):
m4_define(rNumberOfCells, 50)
m4_define(x1NumberOfCells, 40)
m4_define(xInnerNumberOfCells, 20)
m4_define(x2NumberOfCells, 200)
m4_define(zNumberOfCells, 1)
m4_define(yOuterNumberOfCells, 30)
m4_define(yInnerNumberOfCells, 20)
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
    (calc(0.5*D) 0 z1)                              vlabel(V0)
    (calc(0.5*D2) 0 z1)                             vlabel(V1)
    (x2 0 z1)                                       vlabel(V2)
    (x2 calc(0.5*D2/sqrt2) z1)                      vlabel(V3)
    (calc(0.5*D2/sqrt2) calc(0.5*D2/sqrt2) z1)      vlabel(V4)
    (calc(0.5*D/sqrt2) calc(0.5*D/sqrt2) z1)        vlabel(V5)
    (x2 y z1)                                       vlabel(V6)
    (calc(0.5*D2/sqrt2) y z1)                       vlabel(V7)
    (0 y z1)                                        vlabel(V8)
    (0 calc(0.5*D2) z1)                             vlabel(V9)
    (0 calc(0.5*D) z1)                              vlabel(V10)
    (calc(-0.5*D) 0 z1)                             vlabel(V11)
    (calc(-0.5*D2) 0 z1)                            vlabel(V12)
    (x1 0 z1)                                       vlabel(V13)
    (x1 calc(0.5*D2/sqrt2) z1)                      vlabel(V14)
    (-calc(0.5*D2/sqrt2) calc(0.5*D2/sqrt2) z1)     vlabel(V15)
    (-calc(0.5*D/sqrt2) calc(0.5*D/sqrt2) z1)       vlabel(V16)
    (x1 y z1)                                       vlabel(V17)
    (-calc(0.5*D2/sqrt2) y z1)                      vlabel(V18)

    (calc(0.5*D) 0 z2)                              vlabel(V19)
    (calc(0.5*D2) 0 z2)                             vlabel(V20)
    (x2 0 z2)                                       vlabel(V21)
    (x2 calc(0.5*D2/sqrt2) z2)                      vlabel(V22)
    (calc(0.5*D2/sqrt2) calc(0.5*D2/sqrt2) z2)      vlabel(V23)
    (calc(0.5*D/sqrt2) calc(0.5*D/sqrt2) z2)        vlabel(V24)
    (x2 y z2)                                       vlabel(V25)
    (calc(0.5*D2/sqrt2) y z2)                       vlabel(V26)
    (0 y z2)                                        vlabel(V27)
    (0 calc(0.5*D2) z2)                             vlabel(V28)
    (0 calc(0.5*D) z2)                              vlabel(V29)
    (calc(-0.5*D) 0 z2)                             vlabel(V30)
    (calc(-0.5*D2) 0 z2)                            vlabel(V31)
    (x1 0 z2)                                       vlabel(V32)
    (x1 calc(0.5*D2/sqrt2) z2)                      vlabel(V33)
    (-calc(0.5*D2/sqrt2) calc(0.5*D2/sqrt2) z2)     vlabel(V34)
    (-calc(0.5*D/sqrt2) calc(0.5*D/sqrt2) z2)       vlabel(V35)
    (x1 y z2)                                       vlabel(V36)
    (-calc(0.5*D2/sqrt2) y z2)                      vlabel(V37)
);

// Defining blocks:
blocks
(
    hex (5 4 9 10 24 23 28 29) (rNumberOfCells xInnerNumberOfCells zNumberOfCells) simpleGrading (1 1 1)
    hex (0 1 4 5 19 20 23 24) (rNumberOfCells yInnerNumberOfCells zNumberOfCells) simpleGrading (1 1 1)
    hex (1 2 3 4 20 21 22 23) (x2NumberOfCells yInnerNumberOfCells zNumberOfCells) simpleGrading (5 1 1)
    hex (4 3 6 7 23 22 25 26) (x2NumberOfCells yOuterNumberOfCells zNumberOfCells) simpleGrading (5 4 1)
    hex (9 4 7 8 28 23 26 27) (xInnerNumberOfCells yOuterNumberOfCells zNumberOfCells) simpleGrading (1 4 1)
    hex (15 16 10 9 34 35 29 28) (rNumberOfCells xInnerNumberOfCells zNumberOfCells) simpleGrading (1 1 1)
    hex (12 11 16 15 31 30 35 34) (rNumberOfCells yInnerNumberOfCells zNumberOfCells) simpleGrading (1 1 1)
    hex (13 12 15 14 32 31 34 33) (x1NumberOfCells yInnerNumberOfCells zNumberOfCells) simpleGrading (0.2 1 1)
    hex (14 15 18 17 33 34 37 36) (x1NumberOfCells yOuterNumberOfCells zNumberOfCells) simpleGrading (0.2 4 1)
    hex (15 9 8 18 34 28 27 37) (xInnerNumberOfCells yOuterNumberOfCells zNumberOfCells) simpleGrading (1 4 1)
);

edges
(
    arc 0 5 (calc(0.5*D*cos(pi/8)) calc(0.5*D*sin(pi/8)) z1)
    arc 5 10 (calc(0.5*D*cos(3*pi/8)) calc(0.5*D*sin(3*pi/8)) z1)
    arc 1 4 (calc(0.5*D2*cos(pi/8)) calc(0.5*D2*sin(pi/8)) z1)
    arc 4 9 (calc(0.5*D2*cos(3*pi/8)) calc(0.5*D2*sin(3*pi/8)) z1)
    arc 11 16 (calc(-0.5*D*cos(pi/8)) calc(0.5*D*sin(pi/8)) z1)
    arc 16 10 (calc(-0.5*D*cos(3*pi/8)) calc(0.5*D*sin(3*pi/8)) z1)
    arc 12 15 (-calc(0.5*D2*cos(pi/8)) calc(0.5*D2*sin(pi/8)) z1)
    arc 15 9 (-calc(0.5*D2*cos(3*pi/8)) calc(0.5*D2*sin(3*pi/8)) z1)

    arc 19 24 (calc(0.5*D*cos(pi/8)) calc(0.5*D*sin(pi/8)) z2)
    arc 24 29 (calc(0.5*D*cos(3*pi/8)) calc(0.5*D*sin(3*pi/8)) z2)
    arc 20 23 (calc(0.5*D2*cos(pi/8)) calc(0.5*D2*sin(pi/8)) z2)
    arc 23 28 (calc(0.5*D2*cos(3*pi/8)) calc(0.5*D2*sin(3*pi/8)) z2)
    arc 30 35 (calc(-0.5*D*cos(pi/8)) calc(0.5*D*sin(pi/8)) z2)
    arc 35 29 (calc(-0.5*D*cos(3*pi/8)) calc(0.5*D*sin(3*pi/8)) z2)
    arc 31 34 (-calc(0.5*D2*cos(pi/8)) calc(0.5*D2*sin(pi/8)) z2)
    arc 34 28 (-calc(0.5*D2*cos(3*pi/8)) calc(0.5*D2*sin(3*pi/8)) z2)
);

// Defining patches:
boundary
(
    inlet
    {
        type patch;
        faces
        (
            (V13 V14 V33 V32)
            (V14 V17 V36 V33)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (V2 V3 V22 V21)
            (V3 V6 V25 V22)
        );
    }

    cylinder
    {
        type wall;
        faces
        (
            (V0 V5 V24 V19)
            (V5 V10 V29 V24)
            (V10 V16 V35 V29)
            (V16 V11 V30 V35)
        );
    }

    topbottom
    {
        type wall;
        faces
        (
            (V17 V18 V37 V36)
            (V18 V8 V27 V37)
            (V8 V7 V26 V27)
            (V7 V6 V25 V26)
        );
    }


    front
    {
        type empty;
        //neighbourPatch  back;
        faces
        (
            (V13 V12 V15 V14)
            (V12 V11 V16 V15)
            (V16 V10 V9 V15)
            (V10 V5 V4 V9)
            (V0 V1 V4 V5)
            (V1 V2 V3 V4)
            (V14 V15 V18 V17)
            (V15 V9 V8 V18)
            (V9 V4 V7 V8)
            (V4 V3 V6 V7)
        );
    }


    back
    {
        type empty;
        //neighbourPatch  front;
        faces
        (
            (V32 V31 V34 V33)
            (V31 V30 V35 V34)
            (V35 V29 V28 V34)
            (V29 V24 V23 V28)
            (V19 V20 V23 V24)
            (V20 V21 V22 V23)
            (V33 V34 V37 V36)
            (V34 V28 V27 V37)
            (V28 V23 V26 V27)
            (V23 V22 V25 V26)
        );
    }

);

mergePatchPairs
(
);

// ************************************************************************* //
