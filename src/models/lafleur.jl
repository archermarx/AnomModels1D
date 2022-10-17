Base.@kwdef struct LafleurModel <: HallThruster.ZeroEquationModel
    K::Float64 = 4 * sqrt(6)
    α::Float64 = 1.0
end

function HallThruster.initialize_anom!(νan, ::LafleurModel, U, params, i)
    HallThruster.initialize_anom!(νan, HallThruster.TwoZoneBohm, U, params. i)
end

function HallThruster.evaluate_anom(model::LafleurModel, U, params, i)
    (;z_cell, Δz_edge, cache, index, config) = params
    (;K, α) = model

    ncells = size(U, 2)

    if i == ncells
        return cache.νan[end-1]
    elseif i == 1
        return cache.νan[2]
    end

    # upwind differences
    ui_L = U[index.ρiui[1], i-1] / U[index.ρi[1], i-1]
    ui_0 = U[index.ρiui[1], i]   / U[index.ρi[1], i]
    ui_R = U[index.ρiui[1], i+1] / U[index.ρi[1], i+1]

    pe_L = cache.ne[i-1] * cache.Tev[i-1]
    pe_0 = cache.ne[i]   * cache.Tev[i]
    pe_R = cache.ne[i+1] * cache.Tev[i+1]

    (dL, d0, dR) = if ui_0 ≥ 0.0
        coeff = inv(Δz_edge[HallThruster.left_edge(i)])
        (-coeff, coeff, 0.0)
    else
        coeff = inv(Δz_edge[HallThruster.right_edge(i)])
        (0.0, -coeff, coeff)
    end

    ∇_uineTe = dL*ui_L*pe_L + d0*ui_0*pe_0 + dR*ui_R*pe_R

    E = -cache.∇ϕ[i]
    B = cache.B[i]
    ne = cache.ne[i]
    mi = config.propellant.m
    vde = E / B
    cs = sqrt(HallThruster.e * cache.Tev[i] / mi)
    ωce = (HallThruster.e * B / HallThruster.me)

    if params.iteration[] < 2
        νan = ωce / 16
    else
        νan = HallThruster.e * abs(∇_uineTe / (K * HallThruster.me * ne * cs * vde))
    end

    return α * νan + (1 - α) * cache.νan[i]
end
