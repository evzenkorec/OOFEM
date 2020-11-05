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

#include "../sm/CrossSections/Corrosion/simplecorrosioncrosssection.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/Corrosion/corrosionmaterialextensioninterface.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"

namespace oofem {
REGISTER_CrossSection(SimpleCorrosionCrossSection);

  //qqq const  double phaseField --- format problem, gradients of phi and c --- wrong format
void 
SimpleCorrosionCrossSection :: computeStressVector(FloatArray &stress, GaussPoint *gp, const FloatArray &strain, const  double phaseField, TimeStep *tStep)
{
    // This function returns the first Piola-Kirchoff stress in vector format and vector of electrical displacement
    // corresponding to a given deformation gradient according to the stress-deformation
    // mode stored in the each gp.

    MaterialMode mode = gp->giveMaterialMode();
    CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
      if ( !corMat ) {
        OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
      }

    
    if ( mode == _3dMat ) {
      corMat->giveCorrosionRealStressVector(stress, gp, strain, phaseField, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
SimpleCorrosionCrossSection :: computePhaseFieldNfactor(double &N_factor, GaussPoint *gp, double phaseField,  double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
      if ( !corMat ) {
        OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
      }

      corMat->givePhaseField_Nfactor(N_factor, gp, phaseField, gradPhaseField, concentration, gradConcentration, tStep);
}

void
SimpleCorrosionCrossSection :: computePhaseFieldBfactor(double &B_factor, GaussPoint *gp, double phaseField,  double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
      if ( !corMat ) {
        OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
      }

      corMat->givePhaseField_Bfactor(B_factor, gp, phaseField, concentration, gradPhaseField, concentration, gradConcentration, tStep);
}

void
SimpleCorrosionCrossSection :: computeConcentrationNfactor(double &N_factor, GaussPoint *gp, double phaseField, double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
      if ( !corMat ) {
        OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
      }

      corMat->giveConcentration_Nfactor(N_factor, gp, phaseField, concentration, gradPhaseField, concentration, gradConcentration, tStep);
}

void
SimpleCorrosionCrossSection :: computeConcentrationBfactor(double &B_factor, GaussPoint *gp, double phaseField, double concentration, double gradPhaseField, double gradConcentration, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
      if ( !corMat ) {
        OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
      }

      corMat->giveConcentration_Bfactor(B_factor, gp, phaseField, concentration, gradPhaseField, concentration, gradConcentration, tStep);
}

  ///??? dodelat dalsi cleny, PhaseField_Bfactor, ...

  


void
SimpleCorrosionCrossSection :: giveConstitutiveMatrix_uu(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
   CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
   if ( !corMat ) {
     OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
   }

   MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        corMat->giveCorrosion3dMaterialStiffnessMatrix_uu(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

void
SimpleCorrosionCrossSection :: giveConstitutiveMatrix_N_phiphi(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
   CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
   if ( !corMat ) {
     OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
   }

   MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        corMat->give3dMaterialStiffnessMatrix_N_phiphi(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

void
SimpleCorrosionCrossSection :: giveConstitutiveMatrix_B_phiphi(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
   CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
   if ( !corMat ) {
     OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
   }

   MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        corMat->giveCorrosion3dMaterialStiffnessMatrix_B_phiphi(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

void
SimpleCorrosionCrossSection :: giveConstitutiveMatrix_N_cc(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
   CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
   if ( !corMat ) {
     OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
   }

   MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        corMat->give3dMaterialStiffnessMatrix_N_cc(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

void
SimpleCorrosionCrossSection :: giveConstitutiveMatrix_B_cc(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
   CorrosionMaterialExtensionInterface *corMat = static_cast< CorrosionMaterialExtensionInterface * >(this->giveMaterialInterface(CorrosionMaterialExtensionInterfaceType, gp) );
   if ( !corMat ) {
     OOFEM_ERROR("Material doesn't implement the required Corrosion Material interface!");
   }

   MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        corMat->giveCorrosion3dMaterialStiffnessMatrix_B_cc(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

IRResultType
SimpleCorrosionCrossSection :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    double thick = 0.0;
    if ( ir->hasField(_IFT_SimpleCorrosionCrossSection_thick) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thick, _IFT_SimpleCorrosionCrossSection_thick);
        propertyDictionary.add(CS_Thickness, thick);
    }

    double width = 0.0;
    if ( ir->hasField(_IFT_SimpleCorrosionCrossSection_width) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, width, _IFT_SimpleCorrosionCrossSection_width);
        propertyDictionary.add(CS_Width, width);
    }

    double area = 0.0;
    if ( ir->hasField(_IFT_SimpleCorrosionCrossSection_area) ) {
        IR_GIVE_FIELD(ir, area, _IFT_SimpleCorrosionCrossSection_area);
    } else {
        area = thick * width;
    }
    propertyDictionary.add(CS_Area, area);


    this->materialNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_SimpleCorrosionCrossSection_MaterialNumber);
  
    return CrossSection :: initializeFrom(ir);
}


void
SimpleCorrosionCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    Material *mat = domain->giveMaterial(materialNumber);
    MaterialStatus *matStat = mat->CreateStatus(& iGP);
    iGP.setMaterialStatus(matStat);
}


bool
SimpleCorrosionCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    if ( this->giveMaterialNumber() ) {
        return this->domain->giveMaterial( this->giveMaterialNumber() )->isCharacteristicMtrxSymmetric(rMode);
    } else {
        return false; // Bet false...
    }
}


Material
*SimpleCorrosionCrossSection :: giveMaterial(IntegrationPoint *ip)
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


double
SimpleCorrosionCrossSection :: give(int aProperty, GaussPoint *gp)
{
    return this->giveMaterial(gp)->give(aProperty, gp);
}


int
SimpleCorrosionCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }
    return this->giveMaterial(ip)->giveIPValue(answer, ip, type, tStep);
}



int
SimpleCorrosionCrossSection :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    Material *mat = this->giveDomain()->giveMaterial(this->materialNumber);
    
    if ( !dynamic_cast< StructuralMaterial * >(mat) ) {
        OOFEM_WARNING( "material %s without structural support", mat->giveClassName() );
        result = 0;
    }

    return result;
}


void
SimpleCorrosionCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    CrossSection :: giveInputRecord(input);

    if ( this->propertyDictionary.includes(CS_Thickness) ) {
        input.setField(this->propertyDictionary.at(CS_Thickness), _IFT_SimpleCorrosionCrossSection_thick);
    }

    if ( this->propertyDictionary.includes(CS_Width) ) {
        input.setField(this->propertyDictionary.at(CS_Width), _IFT_SimpleCorrosionCrossSection_width);
    }

    if ( this->propertyDictionary.includes(CS_Area) ) {
        input.setField(this->propertyDictionary.at(CS_Area), _IFT_SimpleCorrosionCrossSection_area);
    }

    
    input.setField(this->materialNumber, _IFT_SimpleCorrosionCrossSection_MaterialNumber);
    
}


  
Interface
*SimpleCorrosionCrossSection :: giveMaterialInterface(InterfaceType t, IntegrationPoint *ip)
{
    return this->giveMaterial(ip)->giveInterface(t);
}



int
SimpleCorrosionCrossSection :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->packUnknowns(buff, tStep, gp);
}

int
SimpleCorrosionCrossSection :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
SimpleCorrosionCrossSection :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveMaterial(gp)->estimatePackSize(buff, gp);
}
} // end namespace oofem
