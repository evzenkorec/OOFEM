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


#include "../sm/Elements/Corrosion/baseevzenelement.h"
#include "../sm/Materials/Corrosion/corrosionmaterialextensioninterface.h"

#include "../sm/Materials/structuralms.h"

#include "material.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "domain.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "unknownnumberingscheme.h"


#include <cstdio>

namespace oofem {
BaseEvzenElement :: BaseEvzenElement(int n, Domain *domain)
{
}

    

void
BaseEvzenElement :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_phi, IntArray &locationArray_c, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_phi, const IntArray &dofIdArray_c )
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_phi.clear();
    locationArray_c.clear();

    NLStructuralElement *el = this->giveStructuralElement();
    int k = 0;
    IntArray nodalArray;
    for(int i = 1; i <= el->giveNumberOfDofManagers(); i++) {
      DofManager *dMan = el->giveDofManager( i );
      int itt = 1;
      for(int j = 1; j <= dofIdArray_u.giveSize( ); j++) {
	if(dMan->hasDofID( (DofIDItem) dofIdArray_u.at( j ) )) {
	  //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
	  locationArray_u.followedBy( k + itt);
	}
	itt++;
      }
      for(int j = 1; j <= dofIdArray_phi.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_phi.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_phi.followedBy( k + itt);
	}
	itt++;
      }
      for(int j = 1; j <= dofIdArray_c.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_c.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_c.followedBy( k + itt);
	}
	itt++;
      }
      k += dMan->giveNumberOfDofs( );
    }

    for ( int i = 1; i <= el->giveNumberOfInternalDofManagers(); i++ ) {
      DofManager *dMan = el->giveInternalDofManager( i );
      int itt = 1;
      for(int j = 1; j <= dofIdArray_u.giveSize( ); j++) {
	if(dMan->hasDofID( (DofIDItem) dofIdArray_u.at( j ) )) {
	  //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
	  locationArray_u.followedBy( k + itt);
	  itt++;
	}

      }
      for(int j = 1; j <= dofIdArray_phi.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_phi.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_phi.followedBy( k + itt);
	  itt++;
	}

      }
      for(int j = 1; j <= dofIdArray_c.giveSize( ); j++) {
	if (dMan->hasDofID( (DofIDItem) dofIdArray_c.at( j ) )) {
	  //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
	  locationArray_c.followedBy( k + itt);
	  itt++;
	}

      }

      k += dMan->giveNumberOfDofs( );
    }
    
    

    
}

void
BaseEvzenElement :: computeInternalForcesInputs(FloatArray &stressAnswer, FloatArray strain, FloatArray c, FloatArray phi, FloatArray c_grad, FloatArray phi_grad, FloatArray &pf_Nanswer, FloatArray &pf_Banswer, FloatArray &c_Nanswer, FloatArray &c_Banswer,  GaussPoint *gp, TimeStep *tStep)
{
    //NLStructuralElement *elem = this->giveStructuralElement();
    SimpleCorrosionCrossSection *cs = this->giveCrossSection();
       
    cs->computeStressVector(stressAnswer, gp, strain, phi, tStep);
    cs->computePhaseFieldNfactor(pf_Nanswer, gp, phi, c, phi_grad, c_grad, tStep);
    cs->computePhaseFieldBfactor(pf_Banswer, gp, phi, c, phi_grad, c_grad, tStep);
    cs->computeConcentrationNfactor(c_Nanswer, gp, phi, c, phi_grad, c_grad, tStep);
    cs->computeConcentrationBfactor(c_Banswer, gp, phi, c, phi_grad, c_grad, tStep);
}

void
BaseEvzenElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray d_u;
    FloatMatrix B;
    IntArray IdMask_u;
    this->giveDofManDofIDMask_u( IdMask_u );
    this->giveStructuralElement()->computeVectorOf(IdMask_u, VM_Total, tStep, d_u);
    this->computeBmatrixAt(gp,B);
    answer.beProductOf(B, d_u);
}

void
BaseEvzenElement :: computePhaseField(FloatArray &answer,GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_phi;
    FloatArray d_phi;
    FloatMatrix N_phi;
    this->giveDofManDofIDMask_phi( IdMask_phi );
    this->giveStructuralElement()->computeVectorOf(IdMask_phi, VM_Total, tStep, d_phi);
    this->computePhaseFieldNmatrixAt(gp, N_phi);  
    answer.beProductOf(N_phi,d_phi);
}


void
BaseEvzenElement :: computeConcentration(FloatArray &answer,GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_c;
    FloatArray d_c;
    FloatMatrix N_c;
    this->giveDofManDofIDMask_c( IdMask_c );
    this->giveStructuralElement()->computeVectorOf(IdMask_c, VM_Total, tStep, d_c);
    this->computeConcentrationNmatrixAt(gp, N_c);  
    answer.beProductOf(N_c,d_c);
}

void
BaseEvzenElement :: computePhaseFieldGradient(FloatArray &answer,GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_phi;
    FloatArray d_phi;
    FloatMatrix B_phi;
    this->giveDofManDofIDMask_phi( IdMask_phi );
    this->giveStructuralElement()->computeVectorOf(IdMask_phi, VM_Total, tStep, d_phi);
    this->computePhaseFieldBmatrixAt(gp, B_phi);  
    answer.beProductOf(B_phi,d_phi);
}


void
BaseEvzenElement :: computeConcentrationGradient(FloatArray &answer,GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_c;
    FloatArray d_c;
    FloatMatrix B_c;
    this->giveDofManDofIDMask_c( IdMask_c );
    this->giveStructuralElement()->computeVectorOf(IdMask_c, VM_Total, tStep, d_c);
    this->computeConcentrationBmatrixAt(gp, B_c);  
    answer.beProductOf(B_c,d_c);
}


  



void
BaseEvzenElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveStructuralElement();
    FloatArray BS, stress, strain, c, phi, pf_N, pf_B, c_N, c_B, c_grad, phi_grad, N_phiN, B_phiB, N_cN, B_cB;
    FloatMatrix B_u, N_phi, B_phi, N_c, B_c;
   
    answer.resize(this->giveNumberOfDofs());
    answer.zero();
    
    FloatArray answer_u(this->giveNumberOfDisplacementDofs());
    answer_u.zero();
    FloatArray answer_phi(this->giveNumberOfPhaseFieldDofs());
    answer_phi.zero();
    FloatArray answer_c(this->giveNumberOfConcentrationDofs());
    answer_c.zero();
    
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      this->computeStrainVector(strain, gp, tStep);
      this->computeConcentration(c, gp, tStep);
      this->computePhaseField(phi, gp, tStep);
      this->computeConcentrationGradient(c_grad, gp, tStep);
      this->computePhaseFieldGradient(phi_grad, gp, tStep);
      //
      this->computeInternalForcesInputs(stress, strain, c, phi, c_grad, phi_grad, pf_N, pf_B, c_N, c_B, gp, tstep);
            
      // definition of the function missing
      double dV  = elem->computeVolumeAround(gp);
      // Compute nodal internal forces at nodes as f_u = \int_V B^T*vP dV
      elem->computeBmatrixAt(gp, B_u);  
      BS.beTProductOf(B_u, sigma);
      answer_u.add(dV, BS);
      // Compute nodal internal forces at nodes as f_\phi = \int B^T* vD dV     
      this->computePhaseFieldNmatrixAt(gp, N_phi);
      this->computePhaseFieldBmatrixAt(gp, B_phi);      
      N_phiN.beTProductOf(N_phi, pf_N);
      B_phiB.beTProductOf(B_phi, pf_B);
      answer_phi.add(dV, N_phiN);
      answer_phi.add(dV, B_phiB);
      // Compute nodal internal forces at nodes as f_\phi = \int B^T* vD dV     
      this->computeConcentrationNmatrixAt(gp, N_phi);
      this->computeConcentrationBmatrixAt(gp, B_phi);      
      N_cN.beTProductOf(N_c, c_N);
      B_cB.beTProductOf(B_c, c_B);
      answer_c.add(dV, N_cN);
      answer_c.add(dV, B_cB);
    }

    answer.assemble(answer_u, locationArray_u);
    answer.assemble(answer_phi, locationArray_phi);
    answer.assemble(answer_c, locationArray_c);

}




void
BaseEvzenElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
 
    FloatArray localForces(this->giveNumberOfDisplacementDofs());
    answer.resize(this->giveNumberOfDofs());
    this->computeLocForceLoadVector(localForces, tStep, mode);
    answer.assemble(localForces, locationArray_u);

}


/************************************************************************/
void
BaseEvzenElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveStructuralElement();
    //@todo check this
    //    elem->computeLocalForceLoadVector(answer, tStep, mode);

    // transform result from global cs to nodal cs. if necessary
    if ( answer.isNotEmpty() ) {
        if ( elem->computeGtoLRotationMatrix(T) ) {
            // first back to global cs from element local
            answer.rotatedWith(T, 't');
        }
    } else {
        answer.resize(this->giveNumberOfDisplacementDofs());
        answer.zero();
    }
}


void
BaseEvzenElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    answer.resize(this->giveNumberOfDofs(), this->giveNumberOfDofs());
    answer.zero();


    NLStructuralElement *elem = this->giveStructuralElement();
    SimpleCorrosionCrossSection *cs = this->giveCrossSection(); //unused -- compiler warning ???

    FloatMatrix B_u, B_phi, B_c, N_phi, N_c, D_uu, D_N_phiphi, D_B_phiphi, D_N_cc, D_B_cc, DuuB, DphiphiN, DphiphiB, DccN, DccB;
    FloatMatrix Kuu, Kphiphi, Kcc;
    
    for ( GaussPoint *gp: *elem->giveIntegrationRule(0) ) {
      this->computePhaseFieldBmatrixAt(gp, B_phi);
      this->computeConcentrationBmatrixAt(gp, B_c);

      this->computePhaseFieldNmatrixAt(gp, N_phi);
      this->computeConcentrationNmatrixAt(gp, N_c);
       
      elem->computeBmatrixAt(gp, B_u); 
      cs->giveConstitutiveMatrix_uu(D_uu, rMode, gp, tStep);
      cs->giveConstitutiveMatrix_N_phiphi(D_N_phiphi, rMode, gp, tStep);
      cs->giveConstitutiveMatrix_B_phiphi(D_B_phiphi, rMode, gp, tStep);
      cs->giveConstitutiveMatrix_N_cc(D_N_cc, rMode, gp, tStep);
      cs->giveConstitutiveMatrix_B_cc(D_B_cc, rMode, gp, tStep);
      
      double dV  = elem->computeVolumeAround(gp);
      DuuB.beProductOf(D_uu, B_u);
      DphiphiN.beProductOf(D_N_phiphi, N_phi);
      DphiphiB.beProductOf(D_B_phiphi, B_phi);
      DccN.beProductOf(D_N_cc, N_phi);
      DccB.beProductOf(D_B_cc, B_phi);
            
      Kuu.plusProductUnsym(B_u, DuuB, dV);
      Kphiphi.plusProductUnsym(N_phi, DphiphiN, dV);
      Kphiphi.plusProductUnsym(B_phi, DphiphiB, dV);
      Kcc.plusProductUnsym(N_c, DccN, dV);
      Kcc.plusProductUnsym(B_c, DccB, dV);
    }
      

    answer.assemble(Kuu, locationArray_u);
    answer.assemble(Kphiphi, locationArray_phi);
    answer.assemble(Kcc, locationArray_c);

}

// new name for the crossection ???
SimpleCorrosionCrossSection*
BaseEvzenElement :: giveCrossSection()
// Returns the crossSection of the receiver.
{
  //NLStructuralElement *elem = this->giveElement();
  return static_cast< SimpleCorrosionCrossSection* >( this->giveStructuralElement()->giveCrossSection() );
}



IRResultType
BaseEvzenElement :: initializeFrom(InputRecord *ir)
{
  // @todo Is this function necessary???

    return IRRT_OK;
}

void
BaseEvzenElement :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    FloatArray stress, strain;
    /*
    // force updating strains & stresses
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            this->computeStrainVector(strain, gp, tStep);
            this->computeStressVector(stress, strain, gp, tStep);
        }
    }
    */
}


void
BaseEvzenElement :: postInitialize()
{
  IntArray IdMask_u, IdMask_phi, IdMask_c;
  this->giveDofManDofIDMask_u( IdMask_u );
  this->giveDofManDofIDMask_phi( IdMask_phi );
  this->giveDofManDofIDMask_c( IdMask_c );
  this->giveLocationArrayOfDofIDs(locationArray_u,locationArray_phi,locationArray_c, EModelDefaultEquationNumbering(), IdMask_u, IdMask_phi, IdMask_c);
  
}




} // end namespace oofem


