function [sun_vec_body, sun_vec_eci] = true_sun_vector(utc_time, q)
% TRUE_SUN_VECTOR Returns true sun direction vector in ECI and body frames
%
% Inputs:
%   utc_time - datetime in UTC
%   q        - Quaternion (4x1) from ECI to body
%
% Outputs:
%   sun_vec_body - Unit sun vector in body frame (3x1)
%   sun_vec_eci  - Unit sun vector in ECI frame (3x1)

    jd = juliandate(utc_time);  % Convert to Julian date

    % Get Sun position vector in ECI frame (km)
    rsun_eci_km = planetEphemeris(jd, 'Earth', 'Sun', '405');  % JPL DE405

    % Normalize to get unit vector
    sun_vec_eci = rsun_eci_km' / norm(rsun_eci_km);

    % Convert to body frame using DCM from quaternion
    A_eci2body = quat2dcm(q');
    sun_vec_body = A_eci2body * sun_vec_eci;
end
