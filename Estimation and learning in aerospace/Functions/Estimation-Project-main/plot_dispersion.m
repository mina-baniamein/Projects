function plot_dispersion(w, mag_all, base_color, label_str)
% Plot all FRF magnitudes with low alpha + mean curve.
    if isempty(mag_all), return; end
    hold on;
    % light lines
    n = size(mag_all,1);
    c_light = [base_color, 0.06];  % RGBA for line transparency (R2018a+)
    for k = 1:n
        plot(w, mag_all(k,:), 'Color', c_light, 'HandleVisibility','off');
    end
    % mean curve
    mu = mean(mag_all,1,'omitnan');
    plot(w, mu, '-', 'Color', base_color, 'LineWidth', 1.8, 'DisplayName', label_str);
end

