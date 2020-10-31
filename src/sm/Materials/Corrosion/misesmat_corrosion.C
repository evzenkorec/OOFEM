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

#include "idmgrad.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "sparsemtrx.h"
#include "Materials/isolinearelasticmaterial.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "datastream.h"
#include "contextioerr.h"
#include "stressvector.h"
#include "strainvector.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(MisesCorrosionMaterial);

MisesCorrosionMaterial :: MisesCorrosionMaterial(int n, Domain *d) : MisesMat(n, d), CorrosionMaterialExtensionInterface(d)
//
// constructor
//
{

}


MisesCorrosionMaterial :: ~MisesCorrosionMaterial()
//
// destructor
//
{ }

IRResultType
MisesCorrosionMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    result = MisesMat :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    //inicializace vstupu corrosion, ...    
}



void
MisesCorrosionMaterial :: giveCorrosionRealStressVector(FloatArray &stress, GaussPoint *gp, const FloatArray &totalStrain, double phaseField, TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
  MisesCorrosionMaterialStatus *status = static_cast< MisesCorrosionMaterialStatus * >( this->giveStatus(gp) );

  this->giveRealStressVector(stress, gp, totalStrain, tStep);
  double h = this->computeDegradationFunction(phaseField);
  stress.times(h + kappa);
  
  // update gp
  status->letTempStrainVectorBe(totalStrain);
  status->letTempStressVectorBe(stress);
  
}

void
MisesCorrosionMaterial :: givePhaseField_Nfactor(double &N_factor, GaussPoint *gp, double phaseField, double concentration, TimeStep *tStep)
//
// returns real stress vector in 3d stress space of receiver according to
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
  MisesCorrosionMaterialStatus *status = static_cast< MisesCorrosionMaterialStatus * >( this->giveStatus(gp) );

  double hPrime = this->giveDegradationFunctionPrime(phaseField);
  
  this->giveRealStressVector(stress, gp, totalStrain, tStep);
  double h = this->computeDegradationFunction(phaseField);
  stress.times(h + kappa);
  
  // update gp
  status->letTempStrainVectorBe(totalStrain);
  status->letTempStressVectorBe(stress);
}
  

  


void
VonMisesCorrosionMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                   MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}


void
MisesCorrosionMaterial :: giveCorrosion3dMaterialStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MisesCorrosionMaterialStatus *status = static_cast< VonMisesCorrosionMaterialStatus * >( this->giveStatus(gp) );
  this->giveStiffnessMatrix(answer, mode, gp, tStep);
  double phaseField = status->giveTempPhaseField();
  double h = this->computeDegradationFunction(phaseField);
  answer.times(h + kappa);
}

void
MisesCorrosionMaterial :: giveCorrosion3dMaterialStiffnessMatrix_N_phiphi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MisesCorrosionMaterialStatus *status = static_cast< VonMisesCorrosionMaterialStatus * >( this->giveStatus(gp) );
  
}


  

 



MaterialStatus *
MisesCorrosionMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new MisesCorrosionMaterialStatus(1, VonMisesCorrosionMaterial :: domain, gp);
}


MisesCorrosionMaterialStatus :: MisesCorrosionMaterialStatus(int n, Domain *d, GaussPoint *g) : IsotropicDamageMaterial1Status(n, d, g)
{ }


VonMisesCorrosionMaterialStatus :: ~VonMisesCorrosionMaterialStatus()
{ }




void
VonMisesCorrosionMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterial1Status :: initTempStatus();
    GradientDamageMaterialStatusExtensionInterface :: initTempStatus();
    this->tempDamage = this->damage;

}



void
VonMisesCorrosionMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    IsotropicDamageMaterial1Status :: updateYourself(tStep);
    GradientDamageMaterialStatusExtensionInterface :: updateYourself(tStep);

}



contextIOResultType
VonMisesCorrosionMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    if ( ( iores = IsotropicDamageMaterial1Status :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
VonMisesCorrosionMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    // read parent class status
    if ( ( iores = IsotropicDamageMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


}     // end namespace oofem
