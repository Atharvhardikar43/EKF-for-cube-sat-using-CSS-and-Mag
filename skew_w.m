function skewed_w = skew_w(w)
% SKEW_W Constructs a quaternion Omega matrix from angular velocity vector
%
% Input:
%   w - Angular velocity (3x1) [w1; w2; w3]
% Output:
%   skewed_w - 4x4 matrix used in quaternion derivative computation

    w1 = w(1);
    w2 = w(2);
    w3 = w(3);

    skewed_w = [  0,  -w1,  -w2,  -w3;
                 w1,    0,   w3,  -w2;
                 w2,  -w3,    0,   w1;
                 w3,   w2,  -w1,    0];
end
