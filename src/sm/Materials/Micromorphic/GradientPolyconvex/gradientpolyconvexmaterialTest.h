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

#ifndef gradientpolyconvexTest_h
#define gradientpolyconvexTest_h


#include "../sm/Materials/Micromorphic/micromorphicmaterialextensioninterface.h"
#include "../sm/Materials/Micromorphic/micromorphicms.h"
#include "../sm/Materials/isolinearelasticmaterial.h"
#include "../sm/Materials/structuralmaterial.h"
#include "cltypes.h"


///@name Input fields for MicromorphLEmat
//@{
#define _IFT_GradientPolyconvexMaterialTest_Name "gradientpolyconvexmaterialtest"
#define _IFT_GradientPolyconvexMaterialTest_eps "eps"
#define _IFT_GradientPolyconvexMaterialTest_alpha "alpha"
#define _IFT_GradientPolyconvexMaterialTest_gamma "gamma"
//@}

namespace oofem {

  class GradientPolyconvexMaterialTestStatus : public MicromorphicMaterialStatus
 {
 public:
   GradientPolyconvexMaterialTestStatus(int n, Domain *d, GaussPoint *g, bool sym);  
   ~GradientPolyconvexMaterialTestStatus(){;}
 };

/**
 * MicromorphicLinearElasticMaterial
 */
class GradientPolyconvexMaterialTest : public IsotropicLinearElasticMaterial, MicromorphicMaterialExtensionInterface
{
protected:
  double eps;
  double alpha;
  double gamma;
  FloatMatrix tC1_0, tC2_0;

public:
    GradientPolyconvexMaterialTest(int n, Domain * d);
    virtual ~GradientPolyconvexMaterialTest();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_GradientPolyconvexMaterialTest_Name; }
    virtual const char *giveClassName() const { return "GradientPolyconvexMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == MicromorphicMaterialExtensionInterfaceType ) {
            return static_cast< MicromorphicMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

virtual void giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep){;}
    virtual void giveFiniteStrainGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep);
    virtual void giveFiniteStrainGeneralizedStressVectors_3d (FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep);
    virtual void giveFiniteStrainGeneralizedStressVectors_PlaneStrain (FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep);
    

    void compute_dC_dF(FloatMatrix &dCdF,const FloatArray &vF);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new GradientPolyconvexMaterialTestStatus(1, domain, gp, true); }
    FloatMatrix &givetC1(TimeStep *tStep);
    FloatMatrix &givetC2(TimeStep *tStep);

                                                                     
};


 
 

} // end namespace oofem
#endif
