export
ruleSVBGaussianMeanPrecisionOutNED,
ruleSVBGaussianMeanPrecisionMEND,
ruleMGaussianMeanPrecisionGED,
ruleMGaussianMeanPrecisionEGD


function ruleSVBGaussianMeanPrecisionOutNED(msg_out::Message{F,Univariate},
                                   msg_mean::Message{ExponentialLinearQuadratic},
                                   dist_prec::ProbabilityDistribution) where F<:Gaussian
    dist_mean = msg_mean.dist
    message_prior = ruleSVBGaussianMeanPrecisionOutVGD(nothing, msg_out,dist_prec)
    dist_prior = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance},message_prior.dist)
    approx_dist = dist_prior*msg_mean.dist

    return Message(GaussianMeanVariance, m=unsafeMean(approx_dist), v=unsafeCov(approx_dist))
end

function ruleSVBGaussianMeanPrecisionMEND(msg_out::Message{ExponentialLinearQuadratic},
                                   msg_mean::Message{F, Univariate},
                                   dist_prec::ProbabilityDistribution) where F<:Gaussian

    dist_out = msg_out.dist
    message_prior = ruleSVBGaussianMeanPrecisionOutVGD(nothing, msg_mean,dist_prec)
    dist_prior = convert(ProbabilityDistribution{Univariate, GaussianMeanVariance},message_prior.dist)
    approx_dist = dist_prior*msg_out.dist

    return Message(GaussianMeanVariance, m=unsafeMean(approx_dist), v=unsafeCov(approx_dist))
end
# function ruleSVBGaussianMeanPrecisionOutNED(msg_out::Message{F,Univariate},
#                                    msg_mean::Message{ExponentialLinearQuadratic},
#                                    dist_prec::ProbabilityDistribution) where F<:Gaussian
#
#     msg_mean_prime = approximateDoubleExp(msg_mean)
#     return ruleSVBGaussianMeanPrecisionOutVGD(nothing,msg_mean,dist_prec)
# end
#
# function ruleSVBGaussianMeanPrecisionMEND(msg_out::Message{ExponentialLinearQuadratic},
#                                    msg_mean::Message{F, Univariate},
#                                    dist_prec::ProbabilityDistribution) where F<:Gaussian
#
#     msg_out_prime = approximateDoubleExp(msg_out)
#     return ruleSVBGaussianMeanPrecisionOutVGD(nothing,msg_out_prime,dist_prec)
# end

# function ruleMGaussianMeanPrecisionGED(
#     msg_out::Message{F, Univariate},
#     msg_mean::Message{ExponentialLinearQuadratic},
#     dist_prec::ProbabilityDistribution) where F<:Gaussian
#
#         a = msg_mean.dist.params[:a]
#         b = msg_mean.dist.params[:b]
#         c = msg_mean.dist.params[:c]
#         d = msg_mean.dist.params[:d]
#         m_out, v_out = unsafeMeanCov(msg_out.dist)
#         W_bar = unsafeMean(dist_prec)
#         g(x) = a*x[2]+b*exp(c*x[2] + d*x[2]^2/2)+(x[1]-m_out)^2/v_out + (x[1]-x[2])^2*W_bar
#         msg_mean_prime = approximateDoubleExp(msg_mean)
#         x0 = [m_out; msg_mean_prime.dist.params[:m]]
#         m,Σ = NewtonMethod(g,x0,10)
#
#     return ProbabilityDistribution(Multivariate, GaussianMeanVariance,m=m,v=Σ)
# end
#
function ruleMGaussianMeanPrecisionEGD(
    msg_out::Message{ExponentialLinearQuadratic},
    msg_mean::Message{F, Univariate},
    dist_prec::ProbabilityDistribution) where F<:Gaussian

    a = msg_out.dist.params[:a]
    b = msg_out.dist.params[:b]
    c = msg_out.dist.params[:c]
    d = msg_out.dist.params[:d]
    m_mean, v_mean = unsafeMeanCov(msg_mean.dist)

    W_bar = unsafeMean(dist_prec)
    g(x) = a*x[1]+b*exp(c*x[1] + d*x[1]^2/2)+(x[2]-m_mean)^2/v_mean + (x[1]-x[2])^2*W_bar
    msg_out_prime = approximateDoubleExp(msg_out)
    x0 = [msg_out_prime.dist.params[:m]; m_mean]
    m,Σ = NewtonMethod(g,x0,1)

    return ProbabilityDistribution(Multivariate, GaussianMeanVariance,m=m,v=Σ)
end

# function ruleMGaussianMeanPrecisionEGD(
#     msg_out::Message{ExponentialLinearQuadratic},
#     msg_mean::Message{F, Univariate},
#     dist_prec::ProbabilityDistribution) where F<:Gaussian
#
#     a = msg_out.dist.params[:a]
#     b = msg_out.dist.params[:b]
#     c = msg_out.dist.params[:c]
#     d = msg_out.dist.params[:d]
#     m_mean, v_mean = unsafeMeanCov(msg_mean.dist)
#
#     W_bar = unsafeMean(dist_prec)
#     g(x) = exp(-0.5*(a*x[1]+b*exp(c*x[1] + d*x[1]^2/2)+(x[2]-m_mean)^2/v_mean + (x[1]-x[2])^2*W_bar))
#     m,Σ = multivariateNormalApproximation(g,[-50.0; -50.0],[10; 10])
#
#     return ProbabilityDistribution(Multivariate, GaussianMeanVariance,m=m,v=Σ+1e-8*diageye(2))
# end

# function ruleMGaussianMeanPrecisionGED(
#     msg_out::Message{F, Univariate},
#     msg_mean::Message{ExponentialLinearQuadratic},
#     dist_prec::ProbabilityDistribution) where F<:Gaussian
#
#     msg_mean_prime = approximateDoubleExp(msg_mean)
#
#     return ruleMGaussianMeanPrecisionGGD(msg_out,msg_mean_prime,dist_prec)
# end
#
#
#
# function ruleMGaussianMeanPrecisionEGD(
#     msg_out::Message{ExponentialLinearQuadratic},
#     msg_mean::Message{F, Univariate},
#     dist_prec::ProbabilityDistribution) where F<:Gaussian
#
#     msg_out_prime = approximateDoubleExp(msg_out)
#
#     return ruleMGaussianMeanPrecisionGGD(msg_out_prime,msg_mean,dist_prec)
# end


# ###Custom inbounds
function collectStructuredVariationalNodeInbounds(node::GaussianMeanPrecision, entry::ScheduleEntry, interface_to_msg_idx::Dict{Interface, Int})
    # Collect inbounds
    inbounds = String[]
    entry_recognition_factor_id = recognitionFactorId(entry.interface.edge)
    local_cluster_ids = localRecognitionFactorization(entry.interface.node)

    recognition_factor_ids = Symbol[] # Keep track of encountered recognition factor ids
    for node_interface in entry.interface.node.interfaces
        inbound_interface = ultimatePartner(node_interface)
        node_interface_recognition_factor_id = recognitionFactorId(node_interface.edge)

        if node_interface == entry.interface
            # Ignore marginal of outbound edge
            if (entry.msg_update_rule == SVBGaussianMeanPrecisionOutNED) || (entry.msg_update_rule == SVBGaussianMeanPrecisionMEND)
                inbound_idx = interface_to_msg_idx[inbound_interface]
                push!(inbounds, "messages[$inbound_idx]")
            else
                push!(inbounds, "nothing")
            end
        elseif (inbound_interface != nothing) && isa(inbound_interface.node, Clamp)
            # Hard-code marginal of constant node in schedule
            push!(inbounds, marginalString(inbound_interface.node))
        elseif node_interface_recognition_factor_id == entry_recognition_factor_id
            # Collect message from previous result
            inbound_idx = interface_to_msg_idx[inbound_interface]
            push!(inbounds, "messages[$inbound_idx]")
        elseif !(node_interface_recognition_factor_id in recognition_factor_ids)
            # Collect marginal from marginal dictionary (if marginal is not already accepted)
            marginal_idx = local_cluster_ids[node_interface_recognition_factor_id]
            push!(inbounds, "marginals[:$marginal_idx]")
        end

        push!(recognition_factor_ids, node_interface_recognition_factor_id)
    end

    return inbounds
end
