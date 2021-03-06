/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2015   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef planestrainmicrodil_h
#define planestrainmicrodil_h

#include "../sm/Elements/PlaneStrain/quad1planestrain.h"
#include "../sm/Elements/Micromorphic/basemicromorphicelement.h"

#define _IFT_PlaneStrainMicrodil_Name "planestrainmicrodil"

namespace oofem {
class FEI2dQuadLin;

class PlaneStrainMicrodil : public Quad1PlaneStrain, public BaseMicromorphicElement
{
protected:
    static FEI2dQuadLin interpolation;

public:
    PlaneStrainMicrodil(int n, Domain * d);
    virtual ~PlaneStrainMicrodil() { }

 protected:


    virtual void computeMicromorphicNMatrixAt(GaussPoint *gp, FloatMatrix &N);
    virtual void computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &B);


    virtual NLStructuralElement *giveElement() { return this; }
 public:
    virtual const char *giveInputRecordName() const { return _IFT_PlaneStrainMicrodil_Name; }
    virtual const char *giveClassName() const { return "PlaneStrainMicrodil"; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_m(IntArray &answer);


    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep){BaseMicromorphicElement :: computeStiffnessMatrix(answer, mode, tStep);}

    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep){BaseMicromorphicElement :: computeMassMatrix(answer, tStep);}


    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord){BaseMicromorphicElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);}


    virtual int giveNumberOfMicromorphicDofs(){return 4;}
    virtual int giveNumberOfDisplacementDofs(){return 8;}
    virtual int giveNumberOfDofs(){return 12;}
    virtual void postInitialize();

};
} // end namespace oofem
#endif // planestrainmicrodil_h
