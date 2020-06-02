#ifndef ABLWALLFUNCTION_H
#define ABLWALLFUNCTION_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/core/FieldBCOps.H"

namespace amr_wind {

class ABLWallFunction
{
public:
    ABLWallFunction(const CFDSim& sim);

    ~ABLWallFunction() = default;

    amrex::Real log_law_height() const { return m_log_law_height; }

    amrex::Real utau() const { return m_utau; }

    const amrex::Array<amrex::Real, AMREX_SPACEDIM>& umean() const
    { return m_umean; }

    void init_log_law_height();

    void update_umean(const FieldPlaneAveraging& pa);

    void computeplanar();

    void computeusingheatflux();

    // const amrex::FArrayBox  instplanar() const {return m_store_xy_vel[0];}

    const amrex::FArrayBox&  instplanar() const {return m_store_xy_vel_temp;}

    amrex::Real meanwindspeed() const {return m_mean_windSpeed;}

private:

    const CFDSim& m_sim;
  
    const amrex::AmrCore& m_mesh;

    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_umean{{0.0, 0.0, 0.0}};

    amrex::Real m_log_law_height{0.0};
    amrex::Real m_utau;

    amrex::Real m_kappa{0.41};
    amrex::Real m_z0{0.1};

    amrex::Vector<amrex::Real> m_gravity{{0.0, 0.0, -9.81}};

    int m_direction{2};
    bool m_use_fch{true};

    amrex::Array<amrex::Real, 2> m_coeff_interp{{0.0, 0.0}};

    int m_ncells_x, m_ncells_y;

    amrex::Real m_mean_windSpeed, m_mean_potTemp;

    // amrex::Vector<amrex::Box> m_bx_z_sample;
    // amrex::Vector<amrex::FArrayBox>  m_store_xy_vel;
    // amrex::Vector<int> m_z_sample_index;

    amrex::Box m_bx_z_sample;
    amrex::FArrayBox  m_store_xy_vel_temp;
    int m_z_sample_index;

    bool m_abl_neutral, m_stable, m_heatflux;
    amrex::Real m_abl_surface_flux;
    amrex::Real m_abl_surface_rate;
    amrex::Real m_obhukhov_length;
    amrex::Real m_psi_m, m_psi_h;
    amrex::Real m_abl_surface_temp;

    const amrex::Real m_stable_beta = 5.0;
    const amrex::Real m_unstable_gamma = 16.0;

    const int maxIter=25;

  

};

class ABLVelWallFunc : public FieldBCIface
{
public:
    ABLVelWallFunc(Field& velocity, const ABLWallFunction& wall_func);

    void operator()(Field& velocity, const FieldState rho_state) override;

private:
    const ABLWallFunction& m_wall_func;
};

}

#endif /* ABLWALLFUNCTION_H */
