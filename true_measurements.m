function [sun_vec, mag_vec] = true_measurements(r_ecef, tai_seconds)
    sun_vec = true_sun_vector(tai_seconds);         % From your earlier code
    mag_vec = true_mag_vector(r_ecef, tai_seconds); % From your earlier code
end
