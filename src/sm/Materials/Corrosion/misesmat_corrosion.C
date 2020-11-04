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

    IR_GIVE_FIELD(ir, kappa, _IFT_MisesCorrosionMaterial_kappa); // numerical parameter for convergence
    IR_GIVE_FIELD(ir, D, _IFT_MisesCorrosionMaterial_D); // diffusion koefficient
    IR_GIVE_FIELD(ir, L0, _IFT_MisesCorrosionMaterial_interfaceL0);
    IR_GIVE_FIELD(ir, l, _IFT_MisesCorrosionMaterial_l); // phase filed length scale
    IR_GIVE_FIELD(ir, gamma, _IFT_MisesCorrosionMaterial_gamma); // interface energy
    IR_GIVE_FIELD(ir, cSolid, _IFT_MisesCorrosionMaterial_cSolid); // average concentration of metal
    IR_GIVE_FIELD(ir, cSat, _IFT_MisesCorrosionMaterial_cSat); // average saturation concentration
    IR_GIVE_FIELD(ir, A, _IFT_MisesCorrosionMaterial_A); //free energy density curvature
    IR_GIVE_FIELD(ir, aStar, _IFT_MisesCorrosionMaterial_aStar); //free enrgy density curvature
}

  w = 4.*pow(2, 1/2)*this->aStar*this->gamma/this->l;

  alpha = 16*pow(this->gamma, 2)/w;

  cSe = 1;

  cLe = cSat/cSolid;

double
MisesCorrosionMaterial :: computeInterfaceKineticsCoefficient(){
  return L0;  
}

double
MisesCorrosionMaterial :: computeDegradationFunction(double phaseField){
  return -2.*pow(phaseField, 3) + 3.*pow(phaseField, 2);  
}

double
MisesCorrosionMaterial :: computeDoubleWellPotential(double phaseField){
  return pow(phaseField, 2)*pow(1. - phaseField, 2);  
}

double
MisesCorrosionMaterial :: computeDerivativeOfDegradationFunction(double phaseField){
  return -6.*pow(phaseField, 2) + 6.*phaseField;  
}

double
MisesCorrosionMaterial :: computeSecondDerivativeOfDegradationFunction(double phaseField){
  return -12.*phaseField + 6.;  
}

double
MisesCorrosionMaterial :: computeDerivativeOfDoubleWellPotential(double phaseField){
  return 2.*phaseField*pow(1. - phaseField, 2) - 2.*pow(phaseField, 2)*(1. - phaseField); 
}

double
MisesCorrosionMaterial :: computeSecondDerivativeOfDoubleWellPotential(double phaseField){
  return 2.pow(1. - phaseField, 2) - 8.*phaseField*(1. - phaseField) + 2*pow(phaseField, 2); 
}

double
MisesCorrosionMaterial :: microstress(double phaseField, double concentration){
  double h = this->computeDegradationFunction(phaseField);
  double derOFh = this->computeDerivativeOfDegradationFunction(phaseField);
  double derOFg = this->computeDerivativeOfDoubleWellPotential(phaseField);
  return -2.*A*(concentration - h*(cSe - cLe) - cLe)*(cSe - cLe)*derOFh + w*derOFg;  
}

MisesCorrosionMaterial :: microstressDerivative(double phaseField, double concentration){
  double h = this->computeDegradationFunction(phaseField);
  double derOFh = this->computeDerivativeOfDegradationFunction(phaseField);
  double secondDerOFh = this->computeSecondDerivativeOfDegradationFunction(phaseField);
  double secondDerOFg = this->computeSecondDerivativeOfDoubleWellPotential(phaseField);
  return -2.*A*(-derOFh*(cSe - cLe))*(cSe - cLe)*derOFh - 2.*A*(concentration - h*(cSe - cLe) - cLe)*(cSe - cLe)*secondDerOFh + w*secondDerOFg;  
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
  stress.times(h + this->kappa);
  
  // update gp
  status->letTempStrainVectorBe(totalStrain);
  status->letTempStressVectorBe(stress);
}

void
MisesCorrosionMaterial :: givePhaseField_Nfactor(double &N_factor, GaussPoint *gp, double phaseField, double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)

{
  double h = this->computeDegradationFunction(phaseField);
  double derOFh = this->computeDerivativeOfDegradationFunction(phaseField);
  double micStr = this->microstress(phaseField, concentration);
  double L = this->computeInterfaceKineticsCoefficient();

  N_factor = 1. + L*micStr;
}

void
MisesCorrosionMaterial :: givePhaseField_Bfactor(double &B_factor, GaussPoint *gp, double phaseField,  double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)
{
  double L = this->computeInterfaceKineticsCoefficient();
  B_factor = L*alpha*gradPhaseField;
}

void
MisesCorrosionMaterial :: giveConcentration_Nfactor(double &N_factor, GaussPoint *gp, double phaseField,  double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)
{
  N_factor = 1.;
}

void
MisesCorrosionMaterial :: giveConcentration_Bfactor(double &B_factor, GaussPoint *gp, double phaseField,  double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)
{
  double derOFh = this->computeDerivativeOfDegradationFunction(phaseField); 
  B_factor = D*(gradConcentration - derOFh*(cSe - cLe)*gradPhaseField);
}
  

  


void
MisesCorrosionMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
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
  MisesCorrosionMaterialStatus *status = static_cast< MisesCorrosionMaterialStatus * >( this->giveStatus(gp) );
  //qqq this->giveStiffnessMatrix(answer, mode, gp, tStep);
  double phaseField = status->giveTempPhaseField();
  double h = this->computeDegradationFunction(phaseField);
  answer.times(h + kappa);
}

void
MisesCorrosionMaterial :: giveCorrosion3dMaterialStiffnessMatrix_N_phiphi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MisesCorrosionMaterialStatus *status = static_cast< MisesCorrosionMaterialStatus * >( this->giveStatus(gp) );
  //qqq missing status giveTempPhaseField, giveTempConcentration
  double phaseField = status->giveTempPhaseField();
  double concentration = status->giveTempConcentration();
  double L = this->computeInterfaceKineticsCoefficient();
  double secDerMicStr = this->microstressDerivative(phaseField, concentration);
  answer.times(1. + L*secDerMicStr);
}

void
MisesCorrosionMaterial :: giveCorrosion3dMaterialStiffnessMatrix_B_phiphi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  double L = this->computeInterfaceKineticsCoefficient();  
  answer.times(L*alpha); 
}

void
MisesCorrosionMaterial :: giveCorrosion3dMaterialStiffnessMatrix_N_cc(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  answer.times(1.);   
}
  
void
MisesCorrosionMaterial :: giveCorrosion3dMaterialStiffnessMatrix_B_cc(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  answer.times(D);  
}
 



MaterialStatus *
MisesCorrosionMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new MisesCorrosionMaterialStatus(1, MisesCorrosionMaterial :: domain, gp);
}


MisesCorrosionMaterialStatus :: MisesCorrosionMaterialStatus(int n, Domain *d, GaussPoint *g) : IsotropicDamageMaterial1Status(n, d, g) //qqq
{ }


MisesCorrosionMaterialStatus :: ~MisesCorrosionMaterialStatus()
{ }




void
MisesCorrosionMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    IsotropicDamageMaterial1Status :: initTempStatus();
    GradientDamageMaterialStatusExtensionInterface :: initTempStatus();
    this->tempDamage = this->damage; //qqq

}



void
MisesCorrosionMaterialStatus :: updateYourself(TimeStep *tStep)
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
MisesCorrosionMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
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
MisesCorrosionMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
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
//qqq time dependance missing

//qqq chrasteristic length dependant on time missing

}     // end namespace oofem
