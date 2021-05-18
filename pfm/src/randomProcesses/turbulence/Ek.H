#ifndef Ek_H
#define Ek_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Computes turbulent field following the spectrum of Pope (2000)
// Equation (6.246)
inline tmp<scalarField> Ek
(
    const scalar Cp,
    const scalar p0,
    const scalar beta,
    const scalar nu,
    const scalar cL,
    const scalar cEta,
    const scalar e,
    const scalar epsilon,
    const scalarField& k
)
{
    // Kolmogorov length scale
    const scalar eta = pow(sqr(nu)*nu/epsilon,0.25);
    // Turbulent length scale (Pope, 2000)
    const scalar L = pow(e,1.5)/epsilon;
    
    Info << "L      = " << L   << nl
         << "eta    = " << eta << nl
         << "min(k) = " << min(k) << nl << endl;
    
    tmp<scalarField> tEk
    (
        Cp
       *pow(epsilon,2.0/3.0)
       *pow
        (
            k*L
           /sqrt(sqr(k*L) + cL),
           5.0/3.0 + p0
        )
       *exp
        (
            -beta*(pow(sqr(sqr(k*eta))+sqr(sqr(cEta)),0.25) - cEta)
        )
       /pow(k+vSmall,5.0/3.0)
    );
    Info << "Pope Spectrum computed!" << endl;
    /*
    scalarField& Ekf = tEk();

    label i;
    forAll(Ekf, i)
    {
        if (k[i] < 2 || k[i] > 10)
        {
            Ekf[i] = 0.0;
        }
    }
    */

    return tEk;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
