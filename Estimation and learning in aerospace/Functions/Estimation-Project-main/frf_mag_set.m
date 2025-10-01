function mag_db_all = frf_mag_set(sys_set, w, iy, iu)
    n = numel(sys_set);
    mag_db_all = nan(n, numel(w));
    for k = 1:n
        [mag_db, ~] = frf_mag_phase(sys_set{k}, w, iy, iu);
        mag_db_all(k, :) = mag_db(:).';
    end
end