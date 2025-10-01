function density_complex(z)
% 2D histogram (density map) of complex numbers z on the complex plane.
% Uses histogram2 over Re(z) vs Im(z) and overlays scatter for outliers.
    z = z(~isnan(z) & ~isinf(z));
    if isempty(z)
        text(0.5,0.5,'(no data)','HorizontalAlignment','center'); axis off; return;
    end
    xr = real(z); xi = imag(z);
    if isempty(xr) || isempty(xi)
        text(0.5,0.5,'(no data)','HorizontalAlignment','center'); axis off; return;
    end
    % Binning ranges (auto from data, but padded a bit)
    pad = 0.05;
    xlim_ = [min(xr) max(xr)]; xpad = range_or_one(xlim_) * pad;
    ylim_ = [min(xi) max(xi)]; ypad = range_or_one(ylim_) * pad;
    edgesX = linspace(xlim_(1)-xpad, xlim_(2)+xpad, 60);
    edgesY = linspace(ylim_(1)-ypad, ylim_(2)+ypad, 60);
    histogram2(xr, xi, edgesX, edgesY, 'DisplayStyle','tile', 'Normalization','pdf');
    colorbar; axis tight;
    hold on;
    % Overlay light scatter for visibility of individual points
    plot(xr, xi, '.', 'MarkerSize', 6, 'Color', [0 0 0 0.25], 'HandleVisibility','off');
    % Imag axis line
    yl = ylim; plot([0 0], yl, 'k:', 'HandleVisibility','off');
end
