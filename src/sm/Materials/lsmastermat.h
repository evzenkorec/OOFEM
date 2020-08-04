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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#ifndef lsmastermat_h
#define lsmastermat_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "gaussintegrationrule.h"

///@name Input fields for LargeStrainMasterMaterial
//@{
#define _IFT_LargeStrainMasterMaterial_Name "lsmastermat"
#define _IFT_LargeStrainMasterMaterial_m "m"
#define _IFT_LargeStrainMasterMaterial_slaveMat "slavemat"
#define _IFT_LargeStrainMasterMaterial_InitialStrain "initstrain"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * Large strain master material.
 * Stress and stiffness are computed from small strain(slaveMat) material model
 * using a strain tensor from the Seth-Hill strain tensors family (depends on parameter m,
 * m = 0 logarithmic strain,m = 1 Green-Lagrange strain ...)
 * then stress and stiffness are transformed 2.PK stress and appropriate stiffness
 */
class LargeStrainMasterMaterial : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;
    /// 'slave' material model number.
    int slaveMat;
    /// Specifies the strain tensor.
    double m;
    /// Initial strain
    FloatArray initialStrain;

public:
    LargeStrainMasterMaterial(int n, Domain *d);
    virtual ~LargeStrainMasterMaterial();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveInputRecordName() const { return _IFT_LargeStrainMasterMaterial_Name; }
    virtual const char *giveClassName() const { return "LargeStrainMasterMaterial"; }

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

   virtual void givePlaneStressStiffMtrx_dPdF(FloatMatrix &answer,
                                               MatResponseMode mmode, GaussPoint *gp,
                                               TimeStep *tStep);

    
    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix & answer,
                                                    MatResponseMode,
                                                    GaussPoint * gp,
                                                    TimeStep * tStep);  


    virtual void giveSpatial3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    void giveMembrane2dStiffMtrx_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode, GaussPoint *gp,
                                          TimeStep *tStep);


    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *, const FloatArray &, TimeStep *)
    { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *, const FloatArray &, TimeStep *);
    
    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);
    void giveFirstPKStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &redvF, TimeStep *tStep);
    void giveFirstPKStressVector_Membrane2d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedF, TimeStep *tStep) override;
    
    
    virtual void giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);
    void giveSecondPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

    /// transformation matrices
    void giveTransformationMatrices(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &logStress, const FloatArray &lam, const FloatMatrix &N);
    void giveTransformationMatrices_huhu(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &logStress, const FloatArray &lam, const FloatMatrix &N);
    void giveTransformationMatrices_PlaneStress(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &SetHillStress, const FloatArray &lam, const FloatMatrix &N);
    void giveTransformationMatrices_PlaneStress2(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &SetHillStress, const FloatArray &lam, const FloatMatrix &N);
    void giveTransformationMatrices_PlaneStress_dSdE(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N);
    void giveTransformationMatrices_PlaneStress2_dSdE(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N);

    void giveTransformationMatrices_PlaneStressAnalyticalM_2(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &Fm, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N);
    
    void giveTransformationMatrices_PlaneStressAnalyticalM_2_dSdE(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &Fm, const FloatMatrix &SethHillStress, const FloatArray &lam, const FloatMatrix &N);


    
    
    void giveDeltaS_Product(FloatMatrix &answer, const FloatMatrix &S);
    void giveDeltaS_Product_PlaneStress(FloatMatrix &answer, const FloatMatrix &S);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};

//=============================================================================


class LargeStrainMasterMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatMatrix Pmatrix, TLmatrix, transformationMatrix;
    int slaveMat;
   /// 'slave' gauss point
    GaussPoint *slaveGp;
    /// 'slave' gauss point
    GaussIntegrationRule *irule;

public:
    LargeStrainMasterMaterialStatus(int n, Domain *d, GaussPoint *g, int s);
    virtual ~LargeStrainMasterMaterialStatus();


    const FloatMatrix &givePmatrix() { return Pmatrix; }
    const FloatMatrix &giveTLmatrix() { return TLmatrix; }
    const FloatMatrix &giveTransformationMatrix() { return transformationMatrix; }

    GaussPoint *giveSlaveGp(){return slaveGp;}

    void setPmatrix(const FloatMatrix &values) { Pmatrix = values; }
    void setTLmatrix(const FloatMatrix &values) { TLmatrix = values; }
    void setTransformationMatrix(const FloatMatrix &values) { transformationMatrix = values; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "LargeStrainMasterMaterialStatus"; }
};
} // end namespace oofem
#endif // misesmat_h
