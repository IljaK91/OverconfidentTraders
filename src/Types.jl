@with_kw struct Pars
    # Primitives for any model
    σ²ₐ    = 1.0       # variance of fundamental θ
    σ²ᵤ    = 1.0       # variance of noise traders
    ca     = 0.01      # Parameter for cost function
    cb     = 2.0       # Parameter for cost function
    aₜ     = 0.0       # Mean of distiribution for a
    εₜ     = 0.0       # Mean of distribution for epsilon

    # Switches
    type::Symbol       = :CES   # Payoff function
end

@with_kw struct BPars
    
    ## Economy
    β   # information precision economy
    σ²  # uncertainty economy

    # sliced by signal
    w_z  # weight on z market
    w_zε # weight on ε in the z part
    w_x  # weight on x market
    w_xε # weight on ε in the x=z part
    
    # sliced by fund/noise
    w_p  # weight on prior in the price
    w_ε  # weight on ε in the price
    w_a  # weight on a in the price

    ## Market parameters 
    βᵢ   # information precision market
    σ²ᵢ  # uncertainty market

    # sliced by signal
    wᵢ_z  # weight on z market
    wᵢ_zε # weight on ε in the z part
    wᵢ_x  # weight on x market
    wᵢ_xε # weight on ε in the x=z part
    
    # sliced by fund/noise
    wᵢ_p  # weight on prior in the price
    wᵢ_ε  # weight on ε in the price
    wᵢ_a  # weight on a in the price

    ## Objective parameters 
    σ²ᵢ_obj  # uncertainty objective observer
    
    # sliced by signal
    wᵢ_obj_z  # weight on z objective observer
    wᵢ_obj_zε # weight on ε in the z part
    #wᵢ_obj_x  # weight on x objective observer
    #wᵢ_obj_xε # weight on ε in the x=z part
    
    # sliced by fund/noise
    wᵢ_obj_p  # weight on prior
    wᵢ_obj_ε  # total weight on ε
    wᵢ_obj_a  # total weight on a 

    ## Individual parameters
    βᵢⱼ    # information precision i
    σ²ᵢⱼ   # uncertainty trader i

    # sliced by signal
    wᵢⱼ_z  # weight on z trader i
    wᵢⱼ_zε  # weight on ε in z trader i
    wᵢⱼ_x  # weight on x trader i
    wᵢⱼ_xε # weight on ε in the x=z part 

    # sliced by fund/noise
    wᵢⱼ_p  # weight on prior trader i
    wᵢⱼ_a  # weight on private and public signal combined trader i
    wᵢⱼ_ε # weight on ε in the x=z part 

    # Additional parameters for convenience
    Ωᵢⱼ_x # scaling factor for private signal
    Ωᵢⱼ_z # scaling factor for public signal 
    Ωᵢⱼ_p # scaling factor for prior
    Ωᵢⱼ_ε # scaling factor for noise

    # From the expected utility with linear utility
    σ²ₐi_tilde
    σ²_i_ε_tilde
end