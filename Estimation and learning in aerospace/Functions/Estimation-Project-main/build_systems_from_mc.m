function [systems, mask] = build_systems_from_mc(mc_result, use_converged_only)
% Build a cell array of ss() systems from an MC result struct.
% Returns systems(:) and the logical mask of runs used.
    N = numel(mc_result.A_matrices);
    if isfield(mc_result, 'convergence_flags') && use_converged_only
        mask = mc_result.convergence_flags(:) ~= 0;
    else
        mask = true(N,1);
    end
    idx = find(mask);
    systems = cell(numel(idx),1);
    for k = 1:numel(idx)
        i = idx(k);
        A_ = mc_result.A_matrices{i};
        B_ = mc_result.B_matrices{i};
        C_ = mc_result.C_matrices{i};
        D_ = mc_result.D_matrices{i};
        systems{k} = ss(A_, B_, C_, D_);
    end
end
