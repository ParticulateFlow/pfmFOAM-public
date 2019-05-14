// Parametrized test case for a mold geometry


//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])


//Geometry
// depth of inlet
m4_define(H, 90)
//distance of inlet from top
m4_define(J, 70)
// x-length
m4_define(X, 250)
// y-length
m4_define(Y, 390)
// z-length
m4_define(Z, 70)


//Grid points (integers!):
m4_define(xNumberOfCells, 125)
m4_define(y0NumberOfCells, 100)
m4_define(y1NumberOfCells, 20)
m4_define(y2NumberOfCells, 30)
m4_define(zNumberOfCells, 35)


/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.4.1                                 |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 0.001;

vertices
(
	
	(0 0 0)						vlabel(V0)
	(X 0 0)						vlabel(V1)
	(0 calc(Y-H) 0)					vlabel(V2)
	(X calc(Y-H) 0)					vlabel(V3)
	(0 calc(Y-J) 0)					vlabel(V4)
	(X calc(Y-J) 0)					vlabel(V5)
	(0 Y 0)						vlabel(V6)
	(X Y 0)						vlabel(V7)

	(0 0 Z)						vlabel(V8)
	(X 0 Z)						vlabel(V9)
	(0 calc(Y-H) Z)					vlabel(V10)
	(X calc(Y-H) Z)					vlabel(V11)
	(0 calc(Y-J) Z)					vlabel(V12)
	(X calc(Y-J) Z)					vlabel(V13)
	(0 Y Z)						vlabel(V14)
	(X Y Z)						vlabel(V15)
);

// Defining blocks:
blocks
(
	hex (0 1 3 2 8 9 11 10) (xNumberOfCells y0NumberOfCells zNumberOfCells) simpleGrading (
												( 
													(0.8 0.8 4)
													(0.2 0.2 0.25)
												) 
												0.25 1)
	hex (2 3 5 4 10 11 13 12) (xNumberOfCells y1NumberOfCells zNumberOfCells) simpleGrading (
												( 
													(0.8 0.8 4)
													(0.2 0.2 0.25)
												)
												1 1)
	hex (4 5 7 6 12 13 15 14) (xNumberOfCells y2NumberOfCells zNumberOfCells) simpleGrading (
												( 
													(0.8 0.8 4)
													(0.2 0.2 0.25)
												)
												( 
													(0.5 0.5 2.5)
													(0.5 0.5 0.4)
												) 
												1
												)
);

edges
(
);

// Defining patches:
boundary
(


    frontAndBack
    {
        type wall;
        faces
        (
		(V0 V1 V3 V2)
		(V2 V3 V5 V4)
		(V4 V5 V7 V6)
		(V8 V9 V11 V10)
		(V10 V11 V13 V12)
		(V12 V13 V15 V14)
        );
    }

    leftWall
    {
        type wall;
        faces
        (
		(V0 V8 V10 V2)
		(V2 V10 V12 V4)
		(V4 V12 V14 V6)
        );
    }

    rightWall
    {
        type wall;
        faces
        (
		(V1 V9 V11 V3)
		(V3 V11 V13 V5)
		(V5 V13 V15 V7)
        );
    }

    top
    {
        type wall;
        faces
        (
		(V6 V7 V15 V14)
        );
    }

    outlet
    {
        type patch;
        faces
        (
		(V0 V1 V9 V8)
        );
    }

);

mergePatchPairs 
(
);

// ************************************************************************* //
