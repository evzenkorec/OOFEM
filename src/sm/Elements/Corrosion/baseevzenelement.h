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
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser Base Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser Base Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Base Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef baseeevzenelement_h //why _h ???
#define baseeevzenelement_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/nlstructuralelement.h"

#include "../sm/CrossSections/Corrosion/simplecorrosioncrosssection.h"

namespace oofem {
/**
 * Base class for corrosion formulation.
 * @author Martin Horak
 */
  class BaseEvzenElement
{
protected:
    IntArray displacementDofsOrdering, phaseFieldDofsOrdering, concentrationDofsOrdering;
    IntArray locationArray_u, locationArray_phi, locationArray_c;
    

public:
    BaseEvzenElement(int n, Domain *domain);
    virtual ~BaseEvzenElement() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

protected:
  
    /// Pure virtual functions
    virtual NLStructuralElement *giveStructuralElement() = 0;
    virtual void computePhaseFieldBmatrixAt(GaussPoint *gp, FloatMatrix &Be) = 0;
    virtual void computePhaseFieldNmatrixAt(GaussPoint *gp, FloatMatrix &Ne) = 0;
    virtual void computeConcentrationBmatrixAt(GaussPoint *gp, FloatMatrix &Be) = 0;
    virtual void computeConcentrationNmatrixAt(GaussPoint *gp, FloatMatrix &Ne) = 0;
    virtual int giveNumberOfConcentrationDofs() = 0;
    virtual int giveNumberOfDisplacementDofs() = 0;
    virtual int giveNumberOfPhaseFieldDofs() = 0;
    virtual int giveNumberOfDofs() = 0;

    virtual void giveDofManDofIDMask_u(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_phi(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_c(IntArray &answer) = 0;
    /// End of pure virtual functions

    /// @return Reference to the associated crossSection of element.
    SimpleCorrosionCrossSection *giveCrossSection(); 

    virtual void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void computeInternalForcesInputs(FloatArray &stressAnswer, FloatArray &pf_Nanswer, FloatArray &pf_Banswer, FloatArray &c_Nanswer, FloatArray &c_Banswer, const FloatArray &strain, double c, double phi, const FloatArray &c_grad, const FloatArray &phi_grad,   GaussPoint *gp, TimeStep *tStep);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    void computePhaseField(double &answer,  GaussPoint *gp, TimeStep *tStep);
    void computeConcentration(double &answer,  GaussPoint *gp, TimeStep *tStep);
    void computePhaseFieldGradient(FloatArray &answer,  GaussPoint *gp, TimeStep *tStep);
    void computeConcentrationGradient(FloatArray &answer,  GaussPoint *gp, TimeStep *tStep);

    void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    void computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);

    virtual IntArray &giveDisplacementDofsOrdering() {return displacementDofsOrdering;}
    virtual IntArray &givePhaseFieldDofsOrdering() {return phaseFieldDofsOrdering;}
    virtual IntArray &giveConcentrationDofsOrdering() {return concentrationDofsOrdering;}
    void giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_phi, IntArray &locationArray_c, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_phi,const IntArray &dofIdArray_c);
    virtual void postInitialize();
    virtual void updateInternalState(TimeStep *tStep);


};
} // end namespace oofem

#endif
