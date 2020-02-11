/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    solidEquilibriumDisplacementFoam

Group
    grpStressAnalysisSolvers

Description
    Steady-state segregated finite-volume solver of linear-elastic,
    small-strain deformation of a solid body, with optional thermal
    diffusion and thermal stresses.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D, also generating the
    stress tensor field sigma.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <chrono>
#include <iostream>
#include <stdio.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state segregated finite-volume solver of linear-elastic,"
        " small-strain deformation of a solid body, with optional thermal"
        " diffusion and thermal stresses"
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"

    printf("=== Create Time call\n");
    #include "createTime.H"

    // Calculate create mesh time
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    printf("=== Create Mesh call\n");
    #include "createMesh.H"
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    printf("=== Create Control call\n");
    #include "createControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    // Use for loop for checking mesh reloading
    const int count = 2;
    for (size_t i = 0; i < count; i++) {
        Info<< "\nCalculating displacement field\n" << i << endl;

        // Try to reset timer before run
//        #include "createTime.H"
//        runTime.setTime(0, 0);

        // Calculate time in loops
        std::chrono::steady_clock::time_point loop_begin = std::chrono::steady_clock::now();

        while (runTime.loop())
        {
            Info<< "Iteration: " << runTime.value() << nl << endl;

            #include "readSteadyStressFoamControls.H"

            // Change accFactor
            if (i == 0) {
                accFac = 2;
            } else {
                accFac = 1;
            }
            // Check what is accFac
            Info << "=== accFac is " << accFac << endl;

            solve
            (
                fvm::laplacian(2*mu + lambda, Dcorr, "laplacian(DD,Dcorr)")
              + fvc::div(sigmaExp + sigmaD)
            );

            D += accFac*Dcorr;

            {
                volTensorField gradDcorr(fvc::grad(Dcorr));

                sigmaExp =
                    (lambda - mu)*gradDcorr + mu*gradDcorr.T()
                  + (lambda*I)*tr(gradDcorr);

                sigmaD += accFac*(mu*twoSymm(gradDcorr) + (lambda*I)*tr(gradDcorr));
            }

            #include "calculateStress.H"
            #include "kineticEnergyLimiter.H"

            runTime.printExecutionTime(Info);
        }

        std::chrono::steady_clock::time_point loop_end = std::chrono::steady_clock::now();
        std::cout << "=== Loop " << i << " time difference = " <<
            std::chrono::duration_cast<std::chrono::milliseconds>(loop_end - loop_begin).count() << "[ms]" << std::endl;

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
