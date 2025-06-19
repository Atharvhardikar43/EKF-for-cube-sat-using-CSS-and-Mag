function [B_body, B_eci] = true_mag_vector(position_ecef_m, utc_time, q)
% TRUE_MAG_VECTOR Computes Earth's magnetic field in body and ECI frames
%
% Inputs:
%   position_ecef_m - 3x1 ECEF position [m]
%   utc_time        - datetime in UTC
%   q               - Quaternion (4x1) to rotate ECI to body frame
%
% Outputs:
%   B_body - Normalized magnetic field vector in body frame (3x1)
%   B_eci  - Normalized magnetic field vector in ECI frame (3x1)

    % Convert ECEF to geodetic coordinates
    lla = ecef2lla(position_ecef_m');
    lat = lla(1);
    lon = lla(2);
    alt = lla(3);

    % Convert UTC to decimal year for magnetic model
    utc_time = datetime('now', 'TimeZone', 'UTC');
    utc_time.TimeZone = '';
    decimal_year = decyear(utc_time);

    % Get magnetic field in NED frame (nT)
    B_ned = wrldmagm(alt, lat, lon, decimal_year, '2025');

    % Convert from NED to ECEF (in Tesla)
    [U, V, W] = ned2ecefv(B_ned(1)*1e-9, B_ned(2)*1e-9, B_ned(3)*1e-9, lat, lon);
    B_ecef = [U; V; W];

    % Convert ECEF to ECI
    B_eci = ecef2eci(utc_time, position_ecef_m', B_ecef');

    % Convert ECI to body frame using DCM from quaternion
    A_eci2body = quat2dcm(q');
    B_body = A_eci2body * B_eci(:);

    % Normalize outputs
    B_body = B_body / norm(B_body);
    B_eci = B_eci / norm(B_eci);
end
