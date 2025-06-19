function [q_est, w_est, P_new] = ekf_update(q_est, w_est, P, J, ...
    sun_meas_body, sun_true_body, mag_meas_body, mag_true_body, ...
    sun_true_eci, mag_true_eci, Q, R, dt)
% EKF_UPDATE  Performs one EKF update step for attitude estimation
%
% Inputs:
%   q_est         - Current estimated quaternion (4x1)
%   w_est         - Current estimated angular velocity (3x1)
%   P             - Current state covariance matrix (7x7)
%   J             - Inertia vector [Jx, Jy, Jz]
%   sun_meas_body - Sun sensor measurement in body frame (3x1)
%   sun_true_body - True sun vector in body frame (3x1)
%   mag_meas_body - Magnetometer measurement in body frame (3x1)
%   mag_true_body - True magnetic field in body frame (3x1)
%   sun_true_eci  - True sun vector in inertial frame (3x1)
%   mag_true_eci  - True magnetic field in inertial frame (3x1)
%   Q             - Process noise covariance matrix (7x7)
%   R             - Measurement noise covariance matrix (6x6)
%   dt            - Time step (s)
%
% Outputs:
%   q_est  - Updated quaternion estimate (4x1)
%   w_est  - Updated angular velocity estimate (3x1)
%   P_new  - Updated state covariance matrix (7x7)

    %% 1. Compute Jacobians
    [F, H] = calc_Jacobians(q_est, w_est, J, mag_true_eci, sun_true_eci);

    %% 2. Time update (prediction)
    [q_est, w_est] = propagate_attitude(q_est, w_est, dt);  % Propagate state
    P = F * P * F' + Q;                                     % Predict covariance

    %% 3. Measurement vectors
    z = [mag_meas_body; sun_meas_body];    % Combined sensor measurements (6x1)
    h = [mag_true_body; sun_true_body];    % Expected measurements from state

    %% 4. EKF measurement update
    a = [q_est; w_est];                    % State vector (7x1)
    y = 1 * (z - H * a);                   % Innovation (6x1)
    S = H * P * H' + R;                    % Innovation covariance (6x6)
    K = P * H' / S;                        % Kalman gain (7x6)

    % Compute state update
    update = K * y;                        % 7x1 update vector
    q_est = q_est + update(1:4);          % Apply quaternion correction
    q_est = q_est / norm(q_est);          % Normalize quaternion
    w_est = w_est + update(5:7);          % Apply angular velocity correction

    % Update state covariance
    P_new = (eye(7) - K * H) * P;
end
