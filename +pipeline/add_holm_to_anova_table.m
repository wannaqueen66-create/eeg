function T = add_holm_to_anova_table(T, pCol)
%ADD_HOLM_TO_ANOVA_TABLE Add Holm-adjusted p-values to an ANOVA/fixed-effects table.
%
% T = pipeline.add_holm_to_anova_table(T)
% T = pipeline.add_holm_to_anova_table(T, pCol)
%
% - Detects a p-value column (default: first match among pValue/pvalue/p/ProbF/PrF)
% - Adds/overwrites column: p_holm
% - Keeps NaNs

if nargin < 2 || strlength(string(pCol))==0
    pCol = "";
end

if isempty(T) || ~istable(T)
    return;
end

vars = string(T.Properties.VariableNames);
if pCol==""
    for n = ["pValue","pvalue","p","ProbF","PrF"]
        if any(vars==n)
            pCol = n;
            break;
        end
    end
end

if pCol=="" || ~any(vars==pCol)
    return;
end

try
    p = double(T.(char(pCol)));
catch
    return;
end

try
    T.p_holm = pipeline.holm_stepdown(p);
catch
end
end
