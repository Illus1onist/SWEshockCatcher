/*---------------------------------------------------------------------------*\
    SWEshockFoam: 2D shallow-water solver on an OpenFOAM mesh.

    Equations (flat bed, Cartesian):
        d_t h     + div(h U)               = 0
        d_t (h U) + div(h U U + g h^2/2 I) = 0

    Axisymmetric 1D option (`radial true` in controlDict), r-weighted:
        d_t (r h)   + d_r (r h u)                 = 0
        d_t (r h u) + d_r (r h u^2 + g/2 r h^2)   = g/2 h^2

    Numerics:
      * HLL (1D) / HLLC (2D) Riemann flux; the tangential velocity is
        carried through the HLLC contact;
      * MUSCL reconstruction with a Barth-Jespersen limiter;
      * forward Euler, CFL time step, dry floor and dry tolerance;
      * solid walls / symmetry use F = (0, p_wall n) directly (mirroring
        the state in a Riemann solver adds spurious momentum);
      * radial mode weights the face fluxes by r_f and adds the geometric
        source g/2 h^2.

    The face loop works on an arbitrary mesh; the 1D cases are strips with
    two empty directions, the 2D cases have one empty direction.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// Barth-Jespersen limiter: scale the reconstructed face value so that it
// does not leave the cell-neighbour range.
static Foam::scalar limiter
(
    const Foam::scalar q, const Foam::scalar qMin,
    const Foam::scalar qMax, const Foam::scalar qFace
)
{
    if (qFace > q)
    {
        return Foam::min(1.0, (qMax - q)/Foam::max(qFace - q, SMALL));
    }
    if (qFace < q)
    {
        return Foam::min(1.0, (qMin - q)/Foam::min(qFace - q, -SMALL));
    }
    return 1.0;
}

// Two-wave HLL flux.  Fd is the mass flux and Ff the momentum flux in the
// +nf direction (owner -> neighbour).
static void hll
(
    const Foam::scalar dL, const Foam::vector& UL,
    const Foam::scalar dR, const Foam::vector& UR,
    const Foam::vector& nf, const Foam::scalar g,
    Foam::scalar& Fd, Foam::vector& Ff
)
{
    const Foam::scalar unL = UL & nf;
    const Foam::scalar unR = UR & nf;
    const Foam::scalar cL = Foam::sqrt(g*dL);
    const Foam::scalar cR = Foam::sqrt(g*dR);
    const Foam::scalar SL = Foam::min(unL - cL, unR - cR);
    const Foam::scalar SR = Foam::max(unL + cL, unR + cR);
    const Foam::scalar FLd = dL*unL;
    const Foam::vector FLf = dL*unL*UL + 0.5*g*dL*dL*nf;
    const Foam::scalar FRd = dR*unR;
    const Foam::vector FRf = dR*unR*UR + 0.5*g*dR*dR*nf;

    if (SL >= 0)
    {
        Fd = FLd; Ff = FLf;
    }
    else if (SR <= 0)
    {
        Fd = FRd; Ff = FRf;
    }
    else
    {
        const Foam::scalar den = SR - SL;
        Fd = (SR*FLd - SL*FRd + SL*SR*(dR - dL))/den;
        Ff = (SR*FLf - SL*FRf + SL*SR*(dR*UR - dL*UL))/den;
    }
}

// Three-wave HLLC flux.  The contact/shear wave is needed in 2D; for the
// 1D SWE the two-wave HLL flux is kept (the project's validated choice).
static void hllc
(
    const Foam::scalar dL, const Foam::vector& UL,
    const Foam::scalar dR, const Foam::vector& UR,
    const Foam::vector& nf, const Foam::scalar g,
    const Foam::scalar dMin,
    Foam::scalar& Fd, Foam::vector& Ff
)
{
    const Foam::scalar unL = UL & nf;
    const Foam::scalar unR = UR & nf;
    const Foam::vector utL = UL - unL*nf;
    const Foam::vector utR = UR - unR*nf;
    const Foam::scalar cL = Foam::sqrt(g*dL);
    const Foam::scalar cR = Foam::sqrt(g*dR);
    const Foam::scalar SL = Foam::min(unL - cL, unR - cR);
    const Foam::scalar SR = Foam::max(unL + cL, unR + cR);
    const Foam::scalar pL = 0.5*g*dL*dL;
    const Foam::scalar pR = 0.5*g*dR*dR;

    const Foam::scalar den = dL*(SL - unL) - dR*(SR - unR);
    const Foam::scalar SM = (Foam::mag(den) > SMALL)
        ? (pR - pL + dL*unL*(SL - unL) - dR*unR*(SR - unR))/den
        : 0.5*(unL + unR);

    const Foam::scalar dLs =
        Foam::max(dL*(SL - unL)/Foam::min(SL - SM, -SMALL), dMin);
    const Foam::scalar dRs =
        Foam::max(dR*(SR - unR)/Foam::max(SR - SM, SMALL), dMin);
    const Foam::vector huLs = dLs*(SM*nf + utL);
    const Foam::vector huRs = dRs*(SM*nf + utR);

    const Foam::scalar FLd = dL*unL;
    const Foam::vector FLf = dL*unL*UL + pL*nf;
    const Foam::scalar FRd = dR*unR;
    const Foam::vector FRf = dR*unR*UR + pR*nf;

    if (SL >= 0)
    {
        Fd = FLd; Ff = FLf;
    }
    else if (SM >= 0)
    {
        Fd = FLd + SL*(dLs - dL);
        Ff = FLf + SL*(huLs - dL*UL);
    }
    else if (SR > 0)
    {
        Fd = FRd + SR*(dRs - dR);
        Ff = FRf + SR*(huRs - dR*UR);
    }
    else
    {
        Fd = FRd; Ff = FRf;
    }
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const scalar g = 9.81;
    const scalar dMin = 1e-10;
    const scalar dryTol =
        runTime.controlDict().lookupOrDefault<scalar>("dryTol", 1e-8);
    const scalar cfl =
        runTime.controlDict().lookupOrDefault<scalar>("cfl", 0.1);

    // Below this depth a face falls back to first order: reconstructing
    // h and u separately at a drying front creates spurious velocities.
    const scalar dryRecon =
        runTime.controlDict().lookupOrDefault<scalar>("dryRecon", 1e-4);

    // `radial true` selects the 1D r-weighted axisymmetric form.
    const bool radial =
        runTime.controlDict().lookupOrDefault<bool>("radial", false);
    const scalar maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", 0.01);

    // MUSCL on/off; the first-order fallback is useful for diagnosing
    // boundary and limiter behaviour.
    const bool muscl =
        runTime.controlDict().lookupOrDefault<bool>("muscl", true);

    // More than one geometric dimension -> HLLC; a 1D strip -> HLL.
    const bool multiD = (mesh.nGeometricD() > 1);

    // Optional override: "hll", "hllc" or "auto" (HLL in 1D, HLLC in 2D).
    const word riemann =
        runTime.controlDict().lookupOrDefault<word>("riemann", "auto");
    const bool useHLLC = (riemann == "hllc") || (riemann == "auto" && multiD);

    if (radial && mesh.nGeometricD() != 1)
    {
        FatalErrorInFunction
            << "The radial switch is implemented for the 1D axisymmetric "
            << "strip only (one geometric dimension)." << exit(FatalError);
    }

    volScalarField depth
    (
        IOobject("h", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );
    volVectorField flux
    (
        IOobject("hU", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE),
        mesh
    );

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();
    const scalarField& magSfInternal = mesh.magSf();
    const vectorField& C = mesh.C();
    const vectorField& Cf = mesh.faceCentres();
    const scalarField& V = mesh.V();
    const cellList& cells = mesh.cells();
    const labelListList& cellCells = mesh.cellCells();

    // `mesh.magSf()` is a surfaceScalarField: converting it to a plain
    // scalarField gives only the internal faces.  The cell update and the
    // CFL loop need the areas of *all* faces, boundary faces included.
    scalarField magSf(mesh.nFaces(), 0.0);
    for (label fi = 0; fi < mesh.nInternalFaces(); ++fi)
    {
        magSf[fi] = magSfInternal[fi];
    }
    forAll(mesh.boundary(), pi)
    {
        const fvPatch& p = mesh.boundary()[pi];
        const scalarField& pmag = p.magSf();
        forAll(p, fi)
        {
            magSf[p.start() + fi] = pmag[fi];
        }
    }

    // Radial weighting: r at cell and face centres (axis at x = 0).
    scalarField rCell(mesh.nCells(), 1.0);
    scalarField rFace(mesh.nFaces(), 1.0);
    if (radial)
    {
        forAll(C, i)
        {
            rCell[i] = max(C[i].x(), SMALL);
        }
        forAll(Cf, fi)
        {
            rFace[fi] = max(Cf[fi].x(), 0.0);
        }
    }

    const dimensionedScalar hFloor("hFloor", depth.dimensions(), dMin);

    Info<< "SWEshockFoam: "
        << (radial ? "radial (r-weighted)" : "Cartesian")
        << ", " << (useHLLC ? "HLLC" : "HLL") << " flux" << endl;

    while (runTime.loop())
    {
        scalarField d = depth.internalField();
        vectorField hu = flux.internalField();

        vectorField U = hu/max(d, dMin);
        forAll(d, i)
        {
            d[i] = max(d[i], dMin);
            if (d[i] < dryTol) U[i] = vector::zero;
        }

        // Make the boundary values consistent with the current state before
        // the gradients are evaluated.  `hSafe` avoids 0/0 in hU/h at dry
        // boundaries.
        depth.primitiveFieldRef() = d;
        flux.primitiveFieldRef() = hu;
        depth.correctBoundaryConditions();
        flux.correctBoundaryConditions();

        const volVectorField gradD = fvc::grad(depth);
        const volScalarField hSafe = max(depth, hFloor);
        const volTensorField gradU = fvc::grad(flux/hSafe);

        // Cell min/max for the limiter.
        scalarField dMax = d, dMinC = d;
        vectorField UMax = U, UMin = U;
        forAll(d, i)
        {
            forAll(cellCells[i], j)
            {
                const label m = cellCells[i][j];
                dMax[i] = max(dMax[i], d[m]);
                dMinC[i] = min(dMinC[i], d[m]);
                for (direction k = 0; k < 3; ++k)
                {
                    UMax[i][k] = max(UMax[i][k], U[m][k]);
                    UMin[i][k] = min(UMin[i][k], U[m][k]);
                }
            }
        }

        scalarField Fd(mesh.nFaces(), 0.0);
        vectorField Ff(mesh.nFaces(), vector::zero);

        // ---------- internal faces ----------
        forAll(own, fi)
        {
            const label o = own[fi];
            const label n = nei[fi];
            const vector nf = Sf[fi]/magSf[fi];
            const vector df = 0.5*(C[o] + C[n]) - C[o];
            const vector dn = 0.5*(C[o] + C[n]) - C[n];

            // Unlimited face values.
            scalar dLo, dRo;
            vector ULo, URo;

            if (muscl && d[o] >= dryRecon && d[n] >= dryRecon)
            {
                dLo = d[o] + (gradD[o] & df);
                dRo = d[n] + (gradD[n] & dn);
                ULo = U[o] + (gradU[o] & df);
                URo = U[n] + (gradU[n] & dn);

                // Barth-Jespersen limiter (component-wise for U).
                const scalar phL = limiter(d[o], dMinC[o], dMax[o], dLo);
                const scalar phR = limiter(d[n], dMinC[n], dMax[n], dRo);
                dLo = d[o] + phL*(dLo - d[o]);
                dRo = d[n] + phR*(dRo - d[n]);

                for (direction k = 0; k < 3; ++k)
                {
                    const scalar pl = limiter(U[o][k], UMin[o][k], UMax[o][k], ULo[k]);
                    const scalar pr = limiter(U[n][k], UMin[n][k], UMax[n][k], URo[k]);
                    ULo[k] = U[o][k] + pl*(ULo[k] - U[o][k]);
                    URo[k] = U[n][k] + pr*(URo[k] - U[n][k]);
                }
            }
            else
            {
                dLo = d[o];
                dRo = d[n];
                ULo = U[o];
                URo = U[n];
            }

            const scalar dL = max(dLo, dMin);
            const scalar dR = max(dRo, dMin);

            // Use HLLC only for genuinely wet 2D faces; HLL is more robust
            // at a wet/dry interface.
            if (useHLLC && dL > 10.0*dryTol && dR > 10.0*dryTol)
            {
                hllc(dL, ULo, dR, URo, nf, g, dMin, Fd[fi], Ff[fi]);
            }
            else
            {
                hll(dL, ULo, dR, URo, nf, g, Fd[fi], Ff[fi]);
            }
        }

        // ---------- boundary faces ----------
        forAll(mesh.boundary(), pi)
        {
            const fvPatch& p = mesh.boundary()[pi];
            if (p.type() == "empty") continue;

            const scalarField pd = depth.boundaryField()[pi].patchInternalField();
            const vectorField pu = flux.boundaryField()[pi].patchInternalField();
            const vectorField pSf = p.Sf();
            const scalarField pmag = p.magSf();
            const word& patchType = p.type();
            const word& fieldType = flux.boundaryField()[pi].type();

            const bool wallLike =
                patchType == "wall"
             || patchType == "symmetry"
             || fieldType == "slip"
             || fieldType == "symmetry"
             || fieldType == "wall";

            forAll(p, fi)
            {
                const label facei = p.start() + fi;
                const vector nf = pSf[fi]/max(pmag[fi], SMALL);

                const scalar dL = max(pd[fi], dMin);
                const vector UL = pu[fi]/max(pd[fi], dMin);

                if (wallLike)
                {
                    // Solid wall / symmetry: no mass flux, hydrostatic
                    // pressure only.  This is the exact wall flux; the
                    // mirror-state Riemann trick would add 2h u^2 + h u c.
                    Fd[facei] = 0.0;
                    Ff[facei] = 0.5*g*dL*dL*nf;
                }
                else
                {
                    // zeroGradient / transmissive: ghost state = interior.
                    hll(dL, UL, dL, UL, nf, g, Fd[facei], Ff[facei]);
                }
            }
        }

        // ---------- update ----------
        const scalar dt = runTime.deltaTValue();

        forAll(d, i)
        {
            scalar sumD = 0.0;
            vector sumF = vector::zero;

            // The face flux is stored in the owner -> neighbour direction.
            // For a neighbour cell the outward normal is -nf, hence the
            // sign.  (Adding both faces with + was the instability that
            // blew up the wet pipe run.)
            forAll(cells[i], j)
            {
                const label facei = cells[i][j];
                const scalar sgn = (own[facei] == i) ? 1.0 : -1.0;
                const scalar w = magSf[facei]*rFace[facei];
                sumD += sgn*Fd[facei]*w;
                sumF += sgn*Ff[facei]*w;
            }

            const scalar r = rCell[i];
            scalar q0 = r*d[i] - dt/V[i]*sumD;
            vector q1 = r*hu[i] - dt/V[i]*sumF;

            if (radial)
            {
                // Geometric source of the r-weighted momentum equation.
                q1.x() += dt*0.5*g*d[i]*d[i];
            }

            q0 = max(q0, dMin*r);
            const scalar dn = max(q0/r, dMin);
            vector un = q1/q0;
            if (dn < dryTol) un = vector::zero;

            d[i] = dn;
            U[i] = un;
        }

        forAll(d, i)
        {
            depth[i] = d[i];
            flux[i] = d[i]*U[i];
        }
        depth.correctBoundaryConditions();
        flux.correctBoundaryConditions();

        // CFL: the smallest cell length V/max(A_face), as in OpenFOAM.
        scalar dtNew = GREAT;
        forAll(d, i)
        {
            scalar maxA = 0.0;
            forAll(cells[i], j)
            {
                maxA = max(maxA, magSf[cells[i][j]]);
            }
            const scalar len = V[i]/max(maxA, SMALL);
            const scalar sp = mag(U[i]) + Foam::sqrt(g*d[i]);
            dtNew = min(dtNew, len/max(sp, SMALL));
        }
        runTime.setDeltaT(min(cfl*dtNew, maxDeltaT));

        Info<< "t=" << runTime.timeName()
            << " dt=" << runTime.deltaTValue()
            << " dmin=" << gMin(d) << " dmax=" << gMax(d) << endl;

        runTime.write();
    }

    Info<< "End" << endl;
    return 0;
}
