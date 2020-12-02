"""
    Constructing all Bayesian Weights from the set of parameters and information
    precision.

    Longer than it should be for here, but does the job!
"""
function CstrBPar(βᵢⱼ, βᵢ, β, par::Pars)::BPars
    @unpack_Pars par
    ## Economy
    σ²      = 1/(σ²ₐ^-1 + β + β*σ²ᵤ^-1)

    # sliced by signal
    w_z      = β*σ²ᵤ^-1*σ²
    w_zε     = √β*σ²ᵤ^-1*σ²
    w_x      = β*σ²
    w_xε     = √β*σ²

    # sliced by fund/noise
    w_p      = σ²ₐ^-1*σ²
    w_a      = w_z + w_x
    w_ε      = w_zε + w_xε

    ## Market
    σ²ᵢ      = 1/(σ²ₐ^-1 + βᵢ + βᵢ*σ²ᵤ^-1)

    # sliced by signal
    wᵢ_z      = βᵢ*σ²ᵤ^-1*σ²ᵢ
    wᵢ_zε     = √βᵢ*σ²ᵤ^-1*σ²ᵢ
    wᵢ_x      = βᵢ*σ²ᵢ
    wᵢ_xε     = √βᵢ*σ²ᵢ

    # sliced by fund/noise
    wᵢ_p      = σ²ₐ^-1*σ²ᵢ
    wᵢ_a      = wᵢ_z + wᵢ_x
    wᵢ_ε      = wᵢ_zε + wᵢ_xε

    ## Objective, no private signal!
    σ²ᵢ_obj   = 1/(σ²ₐ^-1 + βᵢ*σ²ᵤ^-1)

    # sliced by signal
    wᵢ_obj_z  = βᵢ*σ²ᵤ^-1*σ²ᵢ_obj
    wᵢ_obj_zε = √βᵢ*σ²ᵤ^-1*σ²ᵢ_obj

    # sliced by fund/noise
    wᵢ_obj_p  = σ²ₐ^-1*σ²ᵢ_obj
    wᵢ_obj_a  = wᵢ_obj_z  
    wᵢ_obj_ε  = wᵢ_obj_zε 

    # Individual
    σ²ᵢⱼ      = 1/(σ²ₐ^-1 + βᵢⱼ + βᵢ*σ²ᵤ^-1)
    # sliced by signal
    wᵢⱼ_z     = βᵢ*σ²ᵤ^-1*σ²ᵢⱼ
    wᵢⱼ_zε    = √βᵢ*σ²ᵤ^-1*σ²ᵢⱼ
    wᵢⱼ_x     = βᵢⱼ*σ²ᵢⱼ
    wᵢⱼ_xε    = √βᵢⱼ*σ²ᵢⱼ

    # sliced by fund/noise
    wᵢⱼ_p     = σ²ₐ^-1*σ²ᵢⱼ
    wᵢⱼ_a     = wᵢⱼ_z  + wᵢⱼ_x
    wᵢⱼ_ε     = wᵢⱼ_zε + wᵢⱼ_xε

    Ωᵢⱼ_x     = wᵢ_x/wᵢⱼ_x
    Ωᵢⱼ_z     = (wᵢ_z - wᵢⱼ_z)/wᵢⱼ_x
    Ωᵢⱼ_p     = (wᵢ_p - wᵢⱼ_p)/wᵢⱼ_x
    Ωᵢⱼ_ε     = (wᵢ_z - wᵢⱼ_z)/wᵢⱼ_xε

    # Perturbed Variances
    σ²ₐi_tilde   = (1/σ²ₐ + βᵢⱼ*(1 - Ωᵢⱼ_z)^2/(1 + σ²ᵤ*Ωᵢⱼ_z^2*βᵢⱼ/β))^-1
    σ²_i_ε_tilde = (1/σ²ᵤ + (Ωᵢⱼ_z^2*βᵢⱼ/β)/(1 + βᵢⱼ*σ²ₐ*(1 - Ωᵢⱼ_z)^2))^-1

    BPars(βᵢ       = βᵢ,
          βᵢⱼ      = βᵢⱼ,
          β        = β,

          ## Economy
          σ²      = σ²,
          w_z     = w_z,
          w_zε    = w_zε,
          w_x     = w_x,
          w_xε    = w_xε,
          w_p     = w_p,
          w_a     = w_a,
          w_ε     = w_ε,

          ### Objective 
          #σ²_obj   = σ²_obj,
          #w_obj_z  = w_obj_z,
          #w_obj_zε = w_obj_zε,

          ## Market
          σ²ᵢ      = σ²ᵢ,
          wᵢ_z     = wᵢ_z,
          wᵢ_zε    = wᵢ_zε,
          wᵢ_x     = wᵢ_x,
          wᵢ_xε    = wᵢ_xε,
          wᵢ_p     = wᵢ_p,
          wᵢ_a     = wᵢ_a,
          wᵢ_ε     = wᵢ_ε,

          ## Objective 
          σ²ᵢ_obj   = σ²ᵢ_obj,
          wᵢ_obj_z  = wᵢ_obj_z,
          wᵢ_obj_zε = wᵢ_obj_zε,

          wᵢ_obj_p  = wᵢ_obj_p,
          wᵢ_obj_a  = wᵢ_obj_a,
          wᵢ_obj_ε  = wᵢ_obj_ε,

          # Individual 
          σ²ᵢⱼ      = σ²ᵢⱼ,
          wᵢⱼ_z     = wᵢⱼ_z,
          wᵢⱼ_zε    = wᵢⱼ_zε,
          wᵢⱼ_x     = wᵢⱼ_x,
          wᵢⱼ_xε    = wᵢⱼ_xε,
          wᵢⱼ_p     = wᵢⱼ_p,
          wᵢⱼ_a     = wᵢⱼ_a,
          wᵢⱼ_ε     = wᵢⱼ_ε,

          Ωᵢⱼ_x     = Ωᵢⱼ_x,
          Ωᵢⱼ_z     = Ωᵢⱼ_z,
          Ωᵢⱼ_p     = Ωᵢⱼ_p,
          Ωᵢⱼ_ε     = Ωᵢⱼ_ε,

          σ²ₐi_tilde   = σ²ₐi_tilde,
          σ²_i_ε_tilde = σ²_i_ε_tilde)
end

# Helper function
CstrBPar(β, par::Pars) = CstrBPar(β, β, β, par)

"""
    Information Acquisition Costs
"""
function infocost(βᵢⱼ,par)
    par.ca*βᵢⱼ^par.cb
end
"""
    Marginal Cost of Information Acquisition
"""
function mc(βᵢⱼ,par)
    par.ca*par.cb*βᵢⱼ^(par.cb-1)
end

# Helper function
mc(βᵢⱼ; par::Pars) = mc(βᵢⱼ,par)

"""
    I only cover the case with an log-linear payoff function. But other cases are also possible!
"""
payoff(a) = exp(a)

"""
    The price function.
"""
function price(z, par, bpar)
    @unpack_Pars par
    @unpack_BPars bpar
    exp(wᵢ_p*aₜ + wᵢ_x*z + wᵢ_z*(z - εₜ/√β) + 1/2*σ²ᵢ)
end

"""
    findbeta(par::Pars; βₑ = 1.)

Finds the symmetric choice of β for which mb(β,β) = mc(β).

# Arguments:
- `par::Pars`: set of parameters
- `βₑ`: keyword argument for market beta for partial equilibrium analysis
"""
function findbeta(par::Pars; βₑ = 1.)

    f = β -> mb(β^2, par) - mc(β^2, par) # Small trick to keep beta positive, but use unbounded rootfinding.
    find_zero(f, 1., atol = 1e-8)^2
end

"""
    Compute the marginal benefit of information acquisition with Gaussian
    Quadrature.
"""
function mb(βᵢⱼ, βᵢ, β, par::Pars; bpar::BPars = CstrBPar(βᵢⱼ, βᵢ, β, par))
    n              = 5
    points, weight = qnwnorm([n,n],[par.aₜ, par.εₜ],[par.σ²ₐ, par.σ²ᵤ])
    a              = points[:,1]
    ϵ              = points[:,2]
    z              = a + ϵ/√βᵢ
    probdiff       = TradeProbDiff.(a, z, par = par, bpar = bpar)
    prof           = profit.(a, z, par = par, bpar = bpar)
    bounds         = 2

    return sum(bounds.*probdiff.*weight.*prof)
end

mb(β, par::Pars) = mb(β, β, β, par)


"""
    Substract the price from the payoff
"""
profit(θ, z; par::Pars, bpar::BPars) = payoff(θ) - price(z, par, bpar)


"""
    Derivative of probability of trading with respect to βᵢ.

    This part depends on the payoff function. Here I provide the analytical expression,
    but in principle it could also be computed numerically. Details can be found in my JMP.
"""
function TradeProbDiff(a, z; par::Pars, bpar::BPars)
    term      = TradeProbTermGrowth(a, z, par, bpar)
    term_diff = TradeProbTermDiffGrowth(a, z, par, bpar)
    term_diff*pdf(Normal(), term)
end

function TradeProbTermGrowth(a, z, par::Pars, bpar::BPars)
    @unpack_Pars par
    @unpack_BPars bpar
    -√βᵢⱼ*(Ωᵢⱼ_p*aₜ + Ωᵢⱼ_x*z + Ωᵢⱼ_z*(z - εₜ/√βᵢ) + 1/(2*wᵢⱼ_x)*(σ²ᵢ - σ²ᵢⱼ) - a)
end

function TradeProbTermDiffGrowth(a, z, par::Pars, bpar::BPars)
    @unpack_Pars par
    @unpack_BPars bpar
    postᵢ          = wᵢ_p*aₜ + wᵢ_x*z + wᵢ_z*(z - εₜ/√βᵢ) + 1/2*σ²ᵢ
    wi_wijx_diff   = -(1/√βᵢⱼ - 1/(2*βᵢⱼ^(3/2))*σ²ᵢⱼ^-1)*postᵢ
    wij_wijx_diff  = -(σ²ₐ^-1*aₜ + βᵢ*σ²ᵤ^-1*(z - εₜ/√βᵢ) + 1/2)/(2*βᵢⱼ^(3/2))

    return term_diff      = 1/(2*√βᵢⱼ)*a + wi_wijx_diff + wij_wijx_diff
end
