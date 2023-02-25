/*!
 * \file CIncNSVariable.hpp
 * \brief Class for defining the variables of the incompressible
          Navier-Stokes solver.
 * \author F. Palacios, T. Economon
 * \version 7.5.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "CIncEulerVariable.hpp"

/*!
 * \class CIncNSVariable
 * \brief Class for defining the variables of the incompressible Navier-Stokes solver.
 * \ingroup Navier_Stokes_Equations
 * \author F. Palacios, T. Economon, T. Albring
 */
class CIncNSVariable final : public CIncEulerVariable {
private:
  VectorType Tau_Wall;        /*!< \brief Magnitude of the wall shear stress from a wall function. */
  VectorType DES_LengthScale;
  su2vector<su2double> Entropy_Generation_Heat_Transfer;
  su2vector<su2double> Entropy_Generation_Viscous_Dissipation;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] pressure - value of the pressure.
   * \param[in] velocity - Value of the flow velocity (initialization value).
   * \param[in] temperature - Value of the temperature (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CIncNSVariable(su2double pressure, const su2double *velocity, su2double temperature,
                 unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config);

  /*!
   * \brief Set the laminar viscosity.
   */
  inline void SetLaminarViscosity(unsigned long iPoint, su2double laminarViscosity) override {
    Primitive(iPoint, indices.LaminarViscosity()) = laminarViscosity;
  }

  /*!
   * \overload
   * \param[in] eddy_visc - Value of the eddy viscosity.
   */
  inline void SetEddyViscosity(unsigned long iPoint, su2double eddy_visc) override {
    Primitive(iPoint, indices.EddyViscosity()) = eddy_visc;
  }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(unsigned long iPoint) const override {
    return Primitive(iPoint, indices.LaminarViscosity());
  }

  /*!
   * \brief Get the eddy viscosity of the flow.
   * \return The eddy viscosity of the flow.
   */
  inline su2double GetEddyViscosity(unsigned long iPoint) const override {
    return Primitive(iPoint, indices.EddyViscosity());
  }

  /*!
   * \brief Set the thermal conductivity.
   */
  inline void SetThermalConductivity(unsigned long iPoint, su2double thermalConductivity) override {
    Primitive(iPoint, indices.ThermalConductivity()) = thermalConductivity;
  }

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(unsigned long iPoint) const override {
    return Primitive(iPoint, indices.ThermalConductivity());
  }

  /*!
   * \brief Set all the primitive variables for incompressible flows
   */
  bool SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel, const su2double *scalar = nullptr);
  using CVariable::SetPrimVar;

  /*!
   * \brief Set the value of the wall shear stress computed by a wall function.
   */
  inline void SetTau_Wall(unsigned long iPoint, su2double tau_wall) override { Tau_Wall(iPoint) = tau_wall; }

  /*!
   * \brief Get the value of the wall shear stress computed by a wall function.
   * \return Value of the wall shear stress computed by a wall function.
   */
  inline su2double GetTau_Wall(unsigned long iPoint) const override { return Tau_Wall(iPoint); }
  inline const VectorType& GetTau_Wall() const { return Tau_Wall; }

  /*!
   * \brief Set the DES Length Scale.
   */
  inline void SetDES_LengthScale(unsigned long iPoint, su2double val_des_lengthscale) override {
    DES_LengthScale(iPoint) = val_des_lengthscale;
  }

  /*!
   * \brief Get the DES length scale
   * \return Value of the DES length Scale.
   */
  inline su2double GetDES_LengthScale(unsigned long iPoint) const override { return DES_LengthScale(iPoint); }

  /*!
   * \brief A virtual member: Set entropy generation due to heat transfer.
   * \param[in] iPoint - Point index.
   * \param[in] val_entropy_gen_ht - value of entropy generation due to heat transfer.
   */
  inline void SetEntropy_Generation_Heat_Transfer(unsigned long iPoint, su2double val_entropy_gen_ht) final {
    Entropy_Generation_Heat_Transfer(iPoint) = val_entropy_gen_ht;
  }

  /*!
   * \brief A virtual member: Get entropy generation due to heat transfer.
   * \param[in] iPoint - Point index.
   * \return Entropy generation due to heat transfer.
   */
  inline su2double GetEntropy_Generation_Heat_Transfer(unsigned long iPoint) const final {
    return Entropy_Generation_Heat_Transfer(iPoint);
  }

  /*!
   * \brief A virtual member: Set entropy generation due to viscous dissipation.
   * \param[in] iPoint - Point index.
   * \param[in] val_entropy_gen_ht - value of entropy generation due to viscous dissipation.
   */
  inline void SetEntropy_Generation_Viscous_Dissipation(unsigned long iPoint, su2double val_entropy_gen_visc) final {
    Entropy_Generation_Viscous_Dissipation(iPoint) = val_entropy_gen_visc;
  }

  /*!
   * \brief A virtual member: Get entropy generation due to viscous dissipation.
   * \param[in] iPoint - Point index.
   * \return Entropy generation due to viscous dissipation.
   */
  inline su2double GetEntropy_Generation_Viscous_Dissipation(unsigned long iPoint) const final {
    return Entropy_Generation_Viscous_Dissipation(iPoint);
  }
};
