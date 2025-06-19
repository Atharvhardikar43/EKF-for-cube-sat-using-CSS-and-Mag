function [F, H] = calc_Jacobians(q_est, w_est, J, b, s)
% CALC_JACOBIANS Computes Jacobians F (dynamics) and H (measurement)
%
% Inputs:
%   q_est - Estimated quaternion (4x1)
%   w_est - Estimated angular velocity (3x1)
%   J     - Inertia vector [Jx, Jy, Jz]
%   b     - Magnetic field vector in inertial frame (3x1)
%   s     - Sun vector in inertial frame (3x1)
%
% Outputs:
%   F     - State transition Jacobian (7x7)
%   H     - Measurement Jacobian (6x7)

    %% 1. Dynamics Jacobian (F)

    W = skew_w(w_est);   % Skew-symmetric matrix from angular velocity
    Q = skew_q(q_est);   % Skew-like matrix from quaternion

    O_3_4 = zeros(3, 4); % 3x4 zero block

    % Angular velocity dynamics matrix (Euler's equation linearized)
    PI = [  0,                          (J(2)-J(3))/J(1)*w_est(3), (J(2)-J(3))/J(1)*w_est(2);
           (J(3)-J(1))/J(2)*w_est(3),  0,                          (J(3)-J(1))/J(2)*w_est(1);
           (J(1)-J(2))/J(3)*w_est(2), (J(1)-J(2))/J(3)*w_est(1),   0 ];

    % State transition Jacobian F: 7x7 matrix
    F = [0.5 * W,     0.5 * Q;
         O_3_4,       PI];

    %% 2. Measurement Jacobian (H)

    % Normalize inertial vectors
    b = b / norm(b);
    s = s / norm(s);

    % Jacobian of magnetic field direction w.r.t quaternion
    B1 = 2 * [ q_est(1),  q_est(4), -q_est(3);
              -q_est(4),  q_est(1),  q_est(2);
               q_est(3), -q_est(2),  q_est(1)] * b;

    B2 = 2 * [ q_est(2),  q_est(3),  q_est(4);
               q_est(3), -q_est(2),  q_est(1);
               q_est(4), -q_est(1), -q_est(2)] * b;

    B3 = 2 * [-q_est(3),  q_est(2), -q_est(1);
               q_est(2),  q_est(3),  q_est(4);
               q_est(1),  q_est(4), -q_est(3)] * b;

    B4 = 2 * [-q_est(4),  q_est(1),  q_est(2);
              -q_est(1), -q_est(4),  q_est(3);
               q_est(2),  q_est(3),  q_est(4)] * b;

    % Jacobian of sun vector direction w.r.t quaternion
    S1 = 2 * [ q_est(1),  q_est(4), -q_est(3);
              -q_est(4),  q_est(1),  q_est(2);
               q_est(3), -q_est(2),  q_est(1)] * s;

    S2 = 2 * [ q_est(2),  q_est(3),  q_est(4);
               q_est(3), -q_est(2),  q_est(1);
               q_est(4), -q_est(1), -q_est(2)] * s;

    S3 = 2 * [-q_est(3),  q_est(2), -q_est(1);
               q_est(2),  q_est(3),  q_est(4);
               q_est(1),  q_est(4), -q_est(3)] * s;

    S4 = 2 * [-q_est(4),  q_est(1),  q_est(2);
              -q_est(1), -q_est(4),  q_est(3);
               q_est(2),  q_est(3),  q_est(4)] * s;

    % Stack all partials into 6x4 matrix
    T = [B1, B2, B3, B4;
         S1, S2, S3, S4];

    O_6_3 = zeros(6, 3); % 6x3 zero block (no measurement dependence on angular velocity)

    % Final measurement Jacobian H: 6x7 matrix
    H = [T, O_6_3];
end
