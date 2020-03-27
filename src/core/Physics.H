#ifndef PHYSICS_H
#define PHYSICS_H

#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"
#include "LevelData.H"

namespace amr_wind {

class Physics
{
public:
    virtual ~Physics() = default;

    virtual std::string physics_name() const = 0;

    //! Initialize the fields for the physics/problem
    virtual void initialize_fields(
        const amrex::Geometry& geom,
        LevelData& leveldata) const = 0;

    //! Perform tasks necessary before advancing timestep
    virtual void pre_advance_work() = 0;

    //! Add momentum source terms
    virtual void add_momentum_sources(
        const amrex::Geometry& geom,
        const LevelData& leveldata,
        amrex::MultiFab& vel_forces) const = 0;

    //! Add momentum source terms
    virtual void add_momentum_sources(
        const amrex::Geometry& geom,
        const amrex::MultiFab& density,
        const amrex::MultiFab& velocity,
        const amrex::MultiFab& scalars_old,
        const amrex::MultiFab& scalars_new,
        amrex::MultiFab& vel_forces) const = 0;
};

}

#endif /* PHYSICS_H */