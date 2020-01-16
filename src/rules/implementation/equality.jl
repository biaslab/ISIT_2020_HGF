export
ruleSPEqualityGaussianGCV

ruleSPEqualityGaussianGCV(msg_1::Message{F1},msg_2::Message{F2},msg_3::Nothing) where {F1<:Gaussian, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_1.dist,msg_2.dist))
ruleSPEqualityGaussianGCV(msg_1::Message{F2},msg_2::Message{F1},msg_3::Nothing) where {F1<:Gaussian, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_1.dist,msg_2.dist))
ruleSPEqualityGaussianGCV(msg_1::Message{F1},msg_2::Nothing,msg_3::Message{F2}) where {F1<:Gaussian, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_1.dist,msg_3.dist))
ruleSPEqualityGaussianGCV(msg_1::Message{F2},msg_2::Nothing,msg_3::Message{F1}) where {F1<:Gaussian, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_1.dist,msg_3.dist))
ruleSPEqualityGaussianGCV(msg_1::Nothing,msg_2::Message{F1},msg_3::Message{F2}) where {F1<:Gaussian, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_2.dist,msg_3.dist))
ruleSPEqualityGaussianGCV(msg_1::Nothing,msg_2::Message{F2},msg_3::Message{F1}) where {F1<:Gaussian, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_2.dist,msg_3.dist))


ruleSPEqualityGaussianGCV(msg_1::Message{F1},msg_2::Message{F2},msg_3::Nothing) where {F1<:ExponentialLinearQuadratic, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_1.dist,msg_2.dist))
ruleSPEqualityGaussianGCV(msg_1::Message{F1},msg_2::Nothing,msg_3::Message{F2}) where {F1<:ExponentialLinearQuadratic, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_1.dist,msg_3.dist))
ruleSPEqualityGaussianGCV(msg_1::Nothing,msg_2::Message{F2},msg_3::Message{F1}) where {F1<:ExponentialLinearQuadratic, F2<:ExponentialLinearQuadratic} = Message(prod!(msg_2.dist,msg_3.dist))
