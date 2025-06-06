module Physical

include("Derivs.jl")
include("Boundary.jl")

#===============================================================================
WaveRHS!:
    * rhs of wave equation
        dot(psi) = Pi
        dot(Pi)  = ddpsi
===============================================================================#
function WaveRHS!(level, r, u)
    psi = u[1]
    Pi = u[2]
    psi_rhs = r[1]
    Pi_rhs = r[2]

    ddpsi = zeros(Float64, level.num_total_points)
    psi_diss = zeros(Float64, level.num_total_points)
    Pi_diss = zeros(Float64, level.num_total_points)
    Derivs.derivs_2nd!(ddpsi, psi, level.dx, level.finite_difference_order)
    Derivs.derivs_diss!(psi_diss, psi, level.dx, level.finite_difference_order)
    Derivs.derivs_diss!(Pi_diss, Pi, level.dx, level.finite_difference_order)

    @. psi_rhs = Pi + level.dissipation * psi_diss
    @. Pi_rhs = ddpsi + level.dissipation * Pi_diss

    if level.is_base_level
        Boundary.ApplyPeriodicBoundaryConditionRHS!(level, r)
    end
end

#===============================================================================
Energy:
    * int_xmin^xmax (Pi^2/2 + dpsi^2/2)
    * calculate on base level (interior) only
===============================================================================#
function Energy(gfs)
    num_total_points = gfs.grid.levels[1].num_total_points
    num_buffer_points = gfs.grid.levels[1].num_buffer_points
    dx = gfs.grid.levels[1].dx
    psi = gfs.levels[1].u[1]
    Pi = gfs.levels[1].u[2]

    dpsi = zeros(Float64, num_total_points)
    Derivs.derivs_1st!(dpsi, psi, dx, gfs.grid.levels[1].finite_difference_order)

    E::Float64 = 0.0
    for i = 1+num_buffer_points:num_total_points-num_buffer_points
        E += (0.5 * Pi[i] * Pi[i] + 0.5 * dpsi[i] * dpsi[i])
    end
    return E * dx
end

end # module Physical
