function [sun_meas_body, sun_meas_eci, mag_meas_body, mag_meas_eci] = ...
    actual_measurements(sun_true_body, sun_true_eci, mag_true_body, mag_true_eci)
% ACTUAL_MEASUREMENTS Simulates sensor readings by adding noise to true vectors
%
% Inputs:
%   sun_true_body - True sun vector in body frame (3x1)
%   sun_true_eci  - True sun vector in ECI frame (3x1)
%   mag_true_body - True magnetic field vector in body frame (3x1)
%   mag_true_eci  - True magnetic field vector in ECI frame (3x1)
%
% Outputs:
%   sun_meas_body - Noisy sun vector in body frame
%   sun_meas_eci  - Noisy sun vector in ECI frame
%   mag_meas_body - Noisy magnetic field in body frame (normalized)
%   mag_meas_eci  - Noisy magnetic field in ECI frame (normalized)

    noise_level = 0.1;  % 1% Gaussian noise

    % Add noise to body-frame vectors
    sun_meas_body = sun_true_body + noise_level * randn(size(sun_true_body));
    mag_meas_body = mag_true_body + noise_level * randn(size(mag_true_body));
    mag_meas_body = mag_meas_body / norm(mag_meas_body);

    % Add noise to inertial-frame vectors
    sun_meas_eci = sun_true_eci + noise_level * randn(size(sun_true_eci));
    mag_meas_eci = mag_true_eci + noise_level * randn(size(mag_true_eci));
    mag_meas_eci = mag_meas_eci / norm(mag_meas_eci);
end
