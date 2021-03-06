#ifndef TaylorGreenVortex_H
#define TaylorGreenVortex_H

#include "Physics.H"
#include "Field.H"

namespace amr_wind {

/** Taylor Green Vortex physics
 */
class TaylorGreenVortex : public Physics::Register<TaylorGreenVortex>
{
public:
    static const std::string identifier() { return "TaylorGreenVortex"; }

    TaylorGreenVortex(const CFDSim& sim);

    virtual ~TaylorGreenVortex() = default;
    //! Initialize the temperature and velocity fields for BoussinesqBubble simulations
    void initialize_fields(
        int level,
        const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:

    Field& m_velocity;
    Field& m_density;

    //! initial density value
    amrex::Real m_rho{1.0};

};

}

#endif /* TaylorGreenVortex_H */
