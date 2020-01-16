@structuredVariationalRule(:node_type => GaussianMeanPrecision,
                           :outbound_type => Message{GaussianMeanVariance},
                           :inbound_types => (Nothing,Message{ExponentialLinearQuadratic},ProbabilityDistribution),
                           :name => SVBGaussianMeanPrecisionOutNED)

@structuredVariationalRule(:node_type => GaussianMeanPrecision,
                          :outbound_type => Message{GaussianMeanVariance},
                          :inbound_types => (Message{ExponentialLinearQuadratic},Nothing,ProbabilityDistribution),
                          :name => SVBGaussianMeanPrecisionMEND)

@marginalRule(:node_type => GaussianMeanPrecision,
            :inbound_types => (Message{Gaussian}, Message{ExponentialLinearQuadratic}, ProbabilityDistribution),
            :name => MGaussianMeanPrecisionGED)

@marginalRule(:node_type => GaussianMeanPrecision,
            :inbound_types => (Message{ExponentialLinearQuadratic}, Message{Gaussian}, ProbabilityDistribution),
            :name => MGaussianMeanPrecisionEGD)
