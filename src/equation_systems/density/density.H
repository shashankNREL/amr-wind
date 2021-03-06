#ifndef DENSITY_H
#define DENSITY_H

#include "PDETraits.H"
#include "SchemeTraits.H"
#include "PDEHelpers.H"
#include "PDE.H"
#include "density/DensitySource.H"

namespace amr_wind {
namespace pde {

struct Density : ScalarTransport
{
    using SrcTerm = DensitySource;

    static std::string pde_name() { return "Density"; }
    static std::string var_name() { return "density"; }

    static constexpr bool multiply_rho = false;
    static constexpr bool has_diffusion = false;

    static constexpr amrex::Real default_bc_value = 1.0;
};

} // namespace pde
} // namespace amr_wind

#endif /* DENSITY_H */
