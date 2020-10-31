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

#ifndef corrosionmaterialextensioninterface_h
#define corrosionmaterialextensioninterface_h

#include "interface.h"
#include "matresponsemode.h"
#include "domain.h"

///@name corrosionmaterialextensioninterface
//@{

//@}

namespace oofem {
class FloatMatrix;
class FloatArray;
class GaussPoint;
class TimeStep;



/**
 * Material interface for gradient material models.
 */
class CorrosionMaterialExtensionInterface : public Interface
{
protected:
    Domain *dom;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    CorrosionMaterialExtensionInterface(Domain *d){    dom = d;}
    /// Destructor.
    virtual ~CorrosionMaterialExtensionInterface() { }

  
    virtual void giveCorrosionRealStressVector_3d(FloatArray &stress, GaussPoint *gp, const FloatArray &strain, double phaseField, TimeStep *tStep) = 0;
  virtual void givePhaseField_Nfactor(double &N_factor, GaussPoint *gp, double phaseField, double concentration, TimeStep *tStep) = 0;
  // givePhaseField_Bfactor
  // ...

  virtual void giveCorrosion3dMaterialStiffnessMatrix_uu(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;
  virtual void giveCorrosion3dMaterialStiffnessMatrix_N_phiphi(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;

virtual void giveCorrosion3dMaterialStiffnessMatrix_B_phiphi(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;

  // virtual void giveCorrosion3dMaterialStiffnessMatrix_N_cc, _B_cc

    
      
};

}
#endif
