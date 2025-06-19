function plot_results(t_log, q_true, q_est, ang_err)
% PLOT_RESULTS Plots quaternion components and angular error over time
%
% Inputs:
%   t_log   - Time vector (datetime array)
%   q_true  - True quaternion log (4 x N)
%   q_est   - Estimated quaternion log (4 x N)
%   ang_err - Angular error (in degrees) over time (1 x N)

    % Compute elapsed time in seconds for plotting
    t_elapsed = seconds(t_log - t_log(1));  

    % Plot each quaternion component
    figure;
    for i = 1:4
        subplot(4, 1, i);
        plot(t_elapsed, q_true(i, :), 'b', ...
             t_elapsed, q_est(i, :), 'r--');
        ylabel(sprintf('q_%d', i - 1));
        grid on;
    end

    % Add legend and labels
    legend('True', 'Estimated');
    xlabel('Time (s)');
    sgtitle('Attitude Estimation: Quaternion Components');

    % Plot angular error
    figure;
    plot(t_elapsed, ang_err);
    xlabel('Time (s)');
    ylabel('Angle Error (deg)');
    title('Attitude Estimation Error');
    grid on;
end
