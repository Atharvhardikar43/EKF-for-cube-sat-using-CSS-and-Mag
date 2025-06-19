function [q_new, w_new] = propagate_attitude(q, w, dt)
% PROPAGATE_ATTITUDE Propagates quaternion forward in time
%
% Inputs:
%   q  - Current quaternion (4x1)
%   w  - Angular velocity vector (3x1) in body frame (rad/s)
%   dt - Time step (s)
%
% Outputs:
%   q_new - Updated quaternion (4x1)
%   w_new - Updated angular velocity (3x1)

    % Construct quaternion kinematic matrix Omega(w)
    Omega = [  0,    -w(1), -w(2), -w(3);
              w(1),   0,     w(3), -w(2);
              w(2), -w(3),   0,     w(1);
              w(3),  w(2), -w(1),   0 ];

    % Quaternion derivative
    q_dot = 0.5 * Omega * q;

    % First-order integration (Euler method)
    q_new = q + q_dot * dt;

    % Normalize to maintain unit quaternion
    q_new = q_new / norm(q_new);

    % No change in angular velocity assumed (can be improved)
    w_new = w;
end
