export HGF

@composite HGF (y, x, u, kappa, omega) begin
    @RV y ~ GaussianMeanVariance(x, exp(kappa*u + omega))
end


@naiveVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                :name          => VariationalHGFOutVPPPP)

@naiveVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                :name          => VariationalHGFIn1PVPPP)

@naiveVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution),
                :name          => VariationalHGFIn2PPVPP)

@naiveVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution),
                :name          => VariationalHGFIn3PPPVP)

@naiveVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing),
                :name          => VariationalHGFIn4PPPPV)

@structuredVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                :name          => SVariationalHGFOutVPPPP)

@structuredVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                :name          => SVariationalHGFIn1PVPPP)

@structuredVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution),
                :name          => SVariationalHGFIn2PPVPP)

@structuredVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution),
                :name          => SVariationalHGFIn3PPPVP)

@structuredVariationalRule(:node_type     => HGF,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing),
                :name          => SVariationalHGFIn4PPPPV)


function ruleVariationalHGFIn1PVPPP(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::Nothing,
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::ProbabilityDistribution{Univariate})


                m_kappa = ForneyLab.unsafeMean(marg_kappa)
                v_kappa = ForneyLab.unsafeCov(marg_kappa)
                m_omega = ForneyLab.unsafeMean(marg_omega)
                v_omega =ForneyLab.unsafeCov(marg_omega)
                m_y = ForneyLab.unsafeMean(marg_y)
                v_y = ForneyLab.unsafeCov(marg_y)
                m_u = ForneyLab.unsafeMean(marg_u)
                v_u = ForneyLab.unsafeCov(marg_u)
                psi = (m_u^2*v_kappa + m_kappa^2*v_u + v_u*v_kappa)
                v_x = v_y + exp((m_u*m_kappa+m_omega - (psi+v_omega)/2))
                Message(GaussianMeanVariance, m=m_y,  v=v_x)

end

function ruleVariationalHGFIn4PPPPV(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::Nothing)


                m_kappa = ForneyLab.unsafeMean(marg_kappa)
                v_kappa = ForneyLab.unsafeCov(marg_kappa)
                m_x = ForneyLab.unsafeMean(marg_x)
                v_x =ForneyLab.unsafeCov(marg_x)
                m_y = ForneyLab.unsafeMean(marg_y)
                v_y = ForneyLab.unsafeCov(marg_y)
                m_u = ForneyLab.unsafeMean(marg_u)
                v_u = ForneyLab.unsafeCov(marg_u)

                a = (m_u^2*v_kappa + m_kappa^2*v_u + v_u*v_kappa)
                psi = ((m_y-m_x)^2+v_y+v_x)*exp(-m_kappa*m_u + a/2)
                # psi_u = clamp(psi, tiny, huge)
                Message(GaussianMeanVariance, m=clamp(log(psi),-huge,huge), v=1.0)
end


function ruleVariationalHGFOutVPPPP(marg_y::Nothing,
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::ProbabilityDistribution{Univariate})


                m_kappa = ForneyLab.unsafeMean(marg_kappa)
                v_kappa = ForneyLab.unsafeCov(marg_kappa)
                m_x = ForneyLab.unsafeMean(marg_x)
                v_x =ForneyLab.unsafeCov(marg_x)
                m_omega = ForneyLab.unsafeMean(marg_omega)
                v_omega = ForneyLab.unsafeCov(marg_omega)
                m_u = ForneyLab.unsafeMean(marg_u)
                v_u = ForneyLab.unsafeCov(marg_u)

                psi = (m_u^2*v_kappa + m_kappa^2*v_u + v_u*v_kappa)
                v_y = v_x + exp((m_u*m_kappa+m_omega - (psi+v_omega)/2))
                Message(GaussianMeanVariance, m=m_x,  v=v_y)
end

function ruleVariationalHGFIn3PPPVP(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::Nothing,
                                    marg_omega::ProbabilityDistribution{Univariate})

                                    m_y = ForneyLab.unsafeMean(marg_y)
                                    v_y = ForneyLab.unsafeCov(marg_y)
                                    m_x = ForneyLab.unsafeMean(marg_x)
                                    v_x =ForneyLab.unsafeCov(marg_x)
                                    m_omega = ForneyLab.unsafeMean(marg_omega)
                                    v_omega = ForneyLab.unsafeCov(marg_omega)
                                    m_u = ForneyLab.unsafeMean(marg_u)
                                    v_u = ForneyLab.unsafeCov(marg_u)
    psi = ((m_y-m_x)^2+v_y+v_x)
    psi_u = clamp(psi, tiny, huge)
    Message(GaussianMeanVariance, m=clamp((log(psi_u)-m_omega+v_omega/2)/m_u,-huge,huge),  v=clamp(2/m_u^2,tiny,huge))
end

function ruleVariationalHGFIn2PPVPP(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::Nothing,
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::ProbabilityDistribution{Univariate})
    m_y = ForneyLab.unsafeMean(marg_y)
    v_y = ForneyLab.unsafeCov(marg_y)
    m_x = ForneyLab.unsafeMean(marg_x)
    v_x =ForneyLab.unsafeCov(marg_x)
    m_omega = ForneyLab.unsafeMean(marg_omega)
    v_omega = ForneyLab.unsafeCov(marg_omega)
    m_kappa = ForneyLab.unsafeMean(marg_kappa)
    v_kappa = ForneyLab.unsafeCov(marg_kappa)
    psi = ((m_y-m_x)^2+v_y+v_x)*exp(-m_omega + v_omega/2)
    psi_u = clamp(psi, tiny, huge)
    Message(GaussianMeanVariance, m=clamp(log(psi_u)/m_kappa , -huge, huge) , v=clamp(2/m_kappa^2, tiny, huge))


end


function ruleSVariationalHGFIn1PVPPP(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::Nothing,
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::ProbabilityDistribution{Univariate})


                m_kappa = ForneyLab.unsafeMean(marg_kappa)
                v_kappa = ForneyLab.unsafeCov(marg_kappa)
                m_omega = ForneyLab.unsafeMean(marg_omega)
                v_omega =ForneyLab.unsafeCov(marg_omega)
                m_y = ForneyLab.unsafeMean(marg_y)
                v_y = ForneyLab.unsafeCov(marg_y)
                m_u = ForneyLab.unsafeMean(marg_u)
                v_u = ForneyLab.unsafeCov(marg_u)
                psi = (m_u^2*v_kappa + m_kappa^2*v_u + v_u*v_kappa)
                v_x = v_y + exp((m_u*m_kappa+m_omega - (psi+v_omega)/2))
                Message(GaussianMeanVariance, m=m_y,  v=v_x)

end

function ruleSVariationalHGFIn4PPPPV(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::Nothing)


                m_kappa = ForneyLab.unsafeMean(marg_kappa)
                v_kappa = ForneyLab.unsafeCov(marg_kappa)
                m_x = ForneyLab.unsafeMean(marg_x)
                v_x =ForneyLab.unsafeCov(marg_x)
                m_y = ForneyLab.unsafeMean(marg_y)
                v_y = ForneyLab.unsafeCov(marg_y)
                m_u = ForneyLab.unsafeMean(marg_u)
                v_u = ForneyLab.unsafeCov(marg_u)

                a = (m_u^2*v_kappa + m_kappa^2*v_u + v_u*v_kappa)
                psi = ((m_y-m_x)^2+v_y+v_x)*exp(-m_kappa*m_u + a/2)
                # psi_u = clamp(psi, tiny, huge)
                Message(GaussianMeanVariance, m=clamp(log(psi),-huge,huge), v=1.0)
end


function ruleSVariationalHGFOutVPPPP(marg_y::Nothing,
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::ProbabilityDistribution{Univariate})


                m_kappa = ForneyLab.unsafeMean(marg_kappa)
                v_kappa = ForneyLab.unsafeCov(marg_kappa)
                m_x = ForneyLab.unsafeMean(marg_x)
                v_x =ForneyLab.unsafeCov(marg_x)
                m_omega = ForneyLab.unsafeMean(marg_omega)
                v_omega = ForneyLab.unsafeCov(marg_omega)
                m_u = ForneyLab.unsafeMean(marg_u)
                v_u = ForneyLab.unsafeCov(marg_u)

                psi = (m_u^2*v_kappa + m_kappa^2*v_u + v_u*v_kappa)
                v_y = v_x + exp((m_u*m_kappa+m_omega - (psi+v_omega)/2))
                Message(GaussianMeanVariance, m=m_x,  v=v_y)
end

function ruleSVariationalHGFIn3PPPVP(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::ProbabilityDistribution{Univariate},
                                    marg_kappa::Nothing,
                                    marg_omega::ProbabilityDistribution{Univariate})

                                    m_y = ForneyLab.unsafeMean(marg_y)
                                    v_y = ForneyLab.unsafeCov(marg_y)
                                    m_x = ForneyLab.unsafeMean(marg_x)
                                    v_x =ForneyLab.unsafeCov(marg_x)
                                    m_omega = ForneyLab.unsafeMean(marg_omega)
                                    v_omega = ForneyLab.unsafeCov(marg_omega)
                                    m_u = ForneyLab.unsafeMean(marg_u)
                                    v_u = ForneyLab.unsafeCov(marg_u)
    psi = ((m_y-m_x)^2+v_y+v_x)
    psi_u = clamp(psi, tiny, huge)
    Message(GaussianMeanVariance, m=clamp((log(psi_u)-m_omega+v_omega/2)/m_u,-huge,huge),  v=clamp(2/m_u^2,tiny,huge))
end

function ruleSVariationalHGFIn2PPVPP(marg_y::ProbabilityDistribution{Univariate},
                                    marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::Nothing,
                                    marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::ProbabilityDistribution{Univariate})
    m_y = ForneyLab.unsafeMean(marg_y)
    v_y = ForneyLab.unsafeCov(marg_y)
    m_x = ForneyLab.unsafeMean(marg_x)
    v_x =ForneyLab.unsafeCov(marg_x)
    m_omega = ForneyLab.unsafeMean(marg_omega)
    v_omega = ForneyLab.unsafeCov(marg_omega)
    m_kappa = ForneyLab.unsafeMean(marg_kappa)
    v_kappa = ForneyLab.unsafeCov(marg_kappa)
    psi = ((m_y-m_x)^2+v_y+v_x)*exp(-m_omega + v_omega/2)
    psi_u = clamp(psi, tiny, huge)
    Message(Gaussian, m=clamp(log(psi_u)/m_kappa , -huge, huge) , v=clamp(2/m_kappa^2, tiny, huge))


end




function ForneyLab.averageEnergy(::Type{HGF}, marg_y::ProbabilityDistribution{Univariate}, marg_x::ProbabilityDistribution{Univariate},
                                    marg_u::ProbabilityDistribution{Univariate}, marg_kappa::ProbabilityDistribution{Univariate},
                                    marg_omega::ProbabilityDistribution{Univariate})

    psi = (unsafeMean(marg_u)^2*unsafeCov(marg_kappa) + unsafeMean(marg_kappa)^2*unsafeCov(marg_u) + unsafeCov(marg_u)*unsafeCov(marg_kappa))
    0.5*log(2*pi) +
    0.5*(unsafeMean(marg_u)*unsafeMean(marg_kappa) + unsafeMean(marg_omega))+

    0.5*(((unsafeMean(marg_y)-unsafeMean(marg_x))^2 + unsafeCov(marg_x) + unsafeCov(marg_y))*exp(-unsafeMean(marg_kappa)*unsafeMean(marg_u)+psi/2)*exp(-unsafeMean(marg_omega)+unsafeCov(marg_omega)/2))

end

# function ForneyLab.averageEnergy(::Type{HGF}, marg_y::ProbabilityDistribution{Univariate}, marg_x::ProbabilityDistribution{Univariate},
#                                     marg_u::ProbabilityDistribution{Univariate}, marg_kappa::ProbabilityDistribution{Univariate,PointMass},
#                                     marg_omega::ProbabilityDistribution{Univariate,PointMass})
#
#     psi =  unsafeCov(marg_u)
#
#     0.5*log(2*pi) +
#     0.5*(unsafeMean(marg_u)+unsafeMean(marg_omega)) +
#     0.5*(((unsafeMean(marg_y)-unsafeMean(marg_x))^2 + unsafeCov(marg_x) + unsafeCov(marg_y))*exp(-unsafeMean(marg_u)+psi/2))*exp(-unsafeMean(marg_omega))
#
# end
