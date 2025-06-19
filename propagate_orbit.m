function [r_new, v_new] = propagate_orbit(r, v, dt, mu, w)
% PROPAGATE_ORBIT Propagates satellite orbit using simple two-body dynamics
%
% Inputs:
%   r   - Current position vector (3x1)
%   v   - Current velocity vector (3x1)
%   dt  - Time step (s)
%   mu  - Gravitational parameter (m^3/s^2)
%   w   - Angular velocity vector (3x1), for rotating frame effects
%
% Outputs:
%   r_new - Updated position (3x1)
%   v_new - Updated velocity (3x1)

    % Gravitational acceleration
    a_grav = -mu * r / norm(r)^3;

    % Coriolis and centrifugal effects (optional, for rotating frame)
    a_coriolis = -2 * cross(w, v);
    a_centrifugal = -cross(w, cross(w, r));

    % Net acceleration (currently using only gravity)
    a = a_grav;

    % Velocity and position update using simple Euler integration
    v_new = v + a * dt;
    r_new = r + v_new * dt;
end
