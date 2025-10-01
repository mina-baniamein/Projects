function [all_poles, all_zeros] = pz_from_systems(sys_set)
% Collect poles and transmission zeros from a set of systems.
    all_poles = [];
    all_zeros = [];
    for k = 1:numel(sys_set)
        sysk = sys_set{k};
        % Poles: eigenvalues of A (same as pole(sys))
        pk = eig(sysk.A);
        all_poles = [all_poles; pk(:)];
        % Transmission zeros: tzero() works for SISO and MIMO
        try
            zk = tzero(sysk);
        catch
            % Fallback: convert to minimal realization first
            zk = tzero(minreal(sysk, 1e-8));
        end
        all_zeros = [all_zeros; zk(:)];
    end
end

