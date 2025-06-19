function skewed_q = skew_q(q)
% SKEW_Q Constructs the quaternion-based linearized matrix used in EKF propagation
%
% Input:
%   q - Quaternion (4x1) [q1; q2; q3; q4]
% Output:
%   skewed_q - 4x3 matrix used in quaternion propagation Jacobian

    q1 = q(1); 
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);

    skewed_q = [-q2, -q3, -q4;
                 q1, -q4,  q3;
                 q4,  q1, -q2;
                -q3,  q2,  q1];
end
