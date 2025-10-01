function [mag_db, ph] = frf_mag_phase(sys, w, iy, iu)
    H = freqresp(sys, w);                % ny × nu × Nw
    Hjw = squeeze(H(iy, iu, :));
    mag_db = 20*log10(abs(Hjw));          % magnitude in dB
    ph     = angle(Hjw);                  % radians
end
