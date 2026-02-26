function export_figure_png(fig, fp, dpi)
%EXPORT_FIGURE_PNG Export figure to PNG without UI toolbars.
%
% Attempts to hide figure/axes toolbars so exported images are clean.
% Uses exportgraphics when available; falls back to saveas.

if nargin < 3 || isempty(dpi)
    dpi = 300;
end

try
    set(fig, 'Color','w');
catch
end

% Hide figure UI
try; fig.ToolBar = 'none'; catch; end
try; fig.MenuBar = 'none'; catch; end

% Hide axes toolbars (newer MATLAB)
try
    axs = findall(fig, 'Type','axes');
    for k=1:numel(axs)
        ax = axs(k);
        try
            axtoolbar(ax, []);
        catch
        end
        try
            if isprop(ax,'Toolbar') && ~isempty(ax.Toolbar)
                ax.Toolbar.Visible = 'off';
            end
        catch
        end
        try
            disableDefaultInteractivity(ax);
        catch
        end
    end
catch
end

% Ensure folder exists
try
    d = fileparts(fp);
    if ~exist(d,'dir'); mkdir(d); end
catch
end

try
    exportgraphics(fig, fp, 'Resolution', dpi);
catch
    try
        saveas(fig, fp);
    catch ME
        fprintf(2,'[WARN] export_figure_png failed: %s\n', ME.message);
    end
end
end
