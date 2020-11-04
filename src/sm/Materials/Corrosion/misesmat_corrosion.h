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

#ifndef electromechnicalms_h
#define electromechnicalms_h

#include "../sm/Materials/structuralms.h"
#include "floatarray.h"
#include "matstatmapperint.h"

#define_IFT_MisesCorrosionMaterial_kappa "kappa"
#define_IFT_MisesCorrosionMaterial_D "diff_coef"
#define_IFT_MisesCorrosionMaterial_interfaceL0 "int_kin_coef_0"
#define_IFT_MisesCorrosionMateria_l "phase_field_length_scale"
#define_IFT_MisesCorrosionMaterial_gamma "interface_energy"
#define_IFT_MisesCorrosionMaterial_cSolid "c_solid"
#define_IFT_MisesCorrosionMaterial_cSat "c_sat"
#define_IFT_MisesCorrosionMaterial_A "free_energy_density_curvature"
#define_IFT_MisesCorrosionMaterial_aStar "astar"


//#define _IFT_VarBasedDamageMaterial_Name "varbaseddamagematerial"

namespace oofem {
class GaussPoint;
class Dictionary;
class Domain;

/**
 * This class implements a structural material status information. It is attribute of
 * gaussPoint. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 *
 * This is a base class for all material statuses corresponding to materials derived from
 * structural material class.
 * It defines stress and strain vectors and their increments.
 * Functions for accessing these components are defined.
 *
 * Tasks:
 * This is abstract class - only basic functionality is supported like:
 * - maintaining and providing access to stress and strain vectors
 *   (including their increments)
 * - storing and restoring status on tape
 * - printingYourself()
 * - updating Yourself after a new equilibrium state has been reached.
 */
class MisesCorrosionMaterial : public StructuralMaterialStatus
{
protected:

  double kappa;
  double D;
  double L0;
  double l;
  double gamma;
  double cSolid;
  double cSat;
  double A;
  double aStar;
  double w;
  double alpha;
  double cSe;
  double cLe;
  

public:
    /// Constructor. Creates new StructuralMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
    MisesCorrosionMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~MisesCorrosionMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual double computeInterfaceKineticsCoefficient();
    virtual double computeDegradationFunction(double phaseField);
    virtual double computeDoubleWellPotential(double phaseField);
    virtual double computeDerivativeOfDegradationFunction(double phaseField);
    virtual double computeDerivativeOfDoubleWellPotential(double phaseField);
    virtual double computeSecondDerivativeOfDegradationFunction(double phaseField);
    virtual double computeSecondDerivativeOfDoubleWellPotential(double phaseField);
    virtual double microstress(double phaseField, double concentration);
    virtual double microstressDerivative(double phaseField, double concentration);

   
    virtual const char *giveClassName() const { return "MisesCorrosionMaterialStatus"; }

    /// Functions for MaterialStatusMapperInterface
    virtual void copyStateVariables(const MaterialStatus &iStatus);
    virtual void addStateVariables(const MaterialStatus &iStatus);
};
} // end namespace oofem
#endif 
