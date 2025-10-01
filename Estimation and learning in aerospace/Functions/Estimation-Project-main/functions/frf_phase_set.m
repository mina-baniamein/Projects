function ph_deg_all = frf_phase_set(sys_set, w, iy, iu)
%FRF_PHASE_SET  Compute phase responses (deg) for a set of systems
%
%   ph_deg_all = frf_phase_set(sys_set, w, iy, iu)
%
% Inputs:
%   sys_set : cell array of LTI systems (ss, tf, etc.)
%   w       : frequency vector [rad/s]
%   iy      : output channel index
%   iu      : input channel index
%
% Output:
%   ph_deg_all : matrix [Nsys x Nfreq] of phase in degrees

    n = numel(sys_set);
    ph_deg_all = nan(n, numel(w));

    for k = 1:n
        sysk = sys_set{k};
        H = freqresp(sysk, w);                  % ny x nu x Nw
        Hjw = squeeze(H(iy, iu, :));            % select channel
        ph_deg = angle(Hjw) * 180/pi;           % convert to degrees
        ph_deg_all(k, :) = ph_deg(:).';
    end
end

