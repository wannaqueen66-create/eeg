function p_adj = holm_stepdown(p)
%HOLM_STEPDOWN Holm step-down p-value adjustment (FWER control).
%
% Usage:
%   p_adj = pipeline.holm_stepdown(p)
%
% Input:
%   p : numeric vector of raw p-values (NaNs allowed)
%
% Output:
%   p_adj : Holm-adjusted p-values (same shape as p)
%
% Notes:
% - Implements Holm (1979) step-down procedure.
% - NaNs are preserved.
% - Values are clamped to [0, 1].

p_adj = nan(size(p));
if isempty(p)
    return;
end

p = double(p);
mask = isfinite(p);
if ~any(mask)
    return;
end

pv = p(mask);
% clamp raw p
pv = max(min(pv,1),0);

m = numel(pv);
[p_sorted, ord] = sort(pv(:), 'ascend');

% raw Holm: p_i * (m - i + 1)
mult = (m:-1:1)';
p_holm = p_sorted .* mult;

% enforce monotonicity (non-decreasing with i)
for i = 2:m
    if p_holm(i) < p_holm(i-1)
        p_holm(i) = p_holm(i-1);
    end
end

% clamp to 1
p_holm = min(p_holm, 1);

% unsort back
p_unsorted = nan(m,1);
p_unsorted(ord) = p_holm;

% write back (linear indexing preserves original shape)
p_adj(mask) = p_unsorted;
end
