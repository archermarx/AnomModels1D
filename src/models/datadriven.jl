Base.@kwdef struct DataDriven1 <: HallThruster.ZeroEquationModel
    c::Vector{Float64} = [2.39, 3.32]
    α::Float64 = 1.0
end

function HallThruster.initialize_anom!(νan, ::DataDriven1, U, params, i)
    HallThruster.initialize_anom!(νan, HallThruster.TwoZoneBohm, U, params. i)
end

function HallThruster.evaluate_anom(model::DataDriven1, U, params, i)
    (;cache, config, index) = params
    (;α, c) = model

    ui = U[index.ρiui[1], i]   / U[index.ρi[1], i]
    E = -cache.∇ϕ[i]
    B = cache.B[i]
    mi = config.propellant.m
    vde = E / B
    cs = sqrt(HallThruster.e * cache.Tev[i] / mi)

    ωce = (HallThruster.e * B / HallThruster.me) / 16

    if params.iteration[] < 2
        νan = ωce / 16
    else
        νan = c[1] * abs(ui / (c[2] * cs + vde)) * ωce
    end

    return α * νan + (1 - α) * cache.νan[i]
end
