function outputDir = plot_nested_M_paper_style(resultRoot, models, outputDir)
%PLOT_NESTED_M_PAPER_STYLE Recreate the paper's 2-by-3 comparison layout.

assert(nargin == 3, 'Use replay_matlab with explicit data and output directories.');
resultRoot = char(resultRoot);
rows = readtable(fullfile(resultRoot, "coverage.csv"), ...
    TextType="string", Delimiter=",");
if ~exist(outputDir, "dir")
    mkdir(outputDir);
end

if nargin < 2
    models = ["linear", "nonlinear", "plate"];
end
methods = ["Student t", "percentile bootstrap", "bootstrap-t"];
colors = [0.10, 0.36, 0.68; 0.80, 0.24, 0.18; 0.10, 0.48, 0.27];
for model = models
    modelRows = rows(rows.model == model, :);
    experimentLabel="independent experiment";bandLabel="95% exact MC interval";
    if model ~= "linear"
        experimentLabel="pool experiment";bandLabel="95% conditional MC interval";
    end
    Mlist = unique(modelRows.M).';
    Rlist = unique(modelRows.R).';
    for M = Mlist
        selected = modelRows(modelRows.M == M, :);
        f = figure("Visible", "off", "Color", "w", ...
            "Units", "centimeters", "Position", [2, 2, 18.2, 10.8]);
        tl = tiledlayout(numel(Rlist), numel(methods), ...
            "TileSpacing", "compact", "Padding", "compact");
        legendHandles = gobjects(3, 1);
        for iR = 1:numel(Rlist)
            R = Rlist(iR);
            rowData = selected(selected.R == R, :);
            rowLimits = localRowLimits(rowData);
            for im = 1:numel(methods)
                method = methods(im);
                ax = nexttile(tl);
                dAll = sortrows(rowData(rowData.method == method, :), "y");
                d = dAll(localDisplayIndices(height(dAll)), :);
                band = fill(ax, [d.y; flipud(d.y)], ...
                    [d.cp95_lower; flipud(d.cp95_upper)], colors(im, :), ...
                    "FaceAlpha", 0.13, "EdgeColor", "none");
                hold(ax, "on");
                formula = plot(ax, dAll.y, dAll.formula_prediction, "-", ...
                    "Color", colors(im, :), "LineWidth", 1.25);
                experiment = plot(ax, d.y, d.coverage, "o", ...
                    "Color", colors(im, :), "MarkerFaceColor", "w", ...
                    "MarkerSize", 3.8, "LineWidth", 0.85);
                yline(ax, 0.95, ":", "Color", [0.32, 0.32, 0.32], ...
                    "LineWidth", 0.75);
                ylim(ax, rowLimits);
                xlim(ax, [min(dAll.y), max(dAll.y)]);
                if iR == 1
                    title(ax, replace(method, "bootstrap-t", "bootstrap t"), "FontWeight", "normal");
                end
                xlabel(ax, localXLabel(model), "Interpreter", "tex");
                if im == 1
                    ylabel(ax, sprintf("R=%d\ncoverage", R));
                else
                    ax.YTickLabel = [];
                end
                box(ax, "on");
                grid(ax, "off");
                ax.FontName = "Times New Roman";
                ax.FontSize = 8.5;
                ax.LineWidth = 0.65;
                ax.TickDir = "in";
                if iR == 1 && im == 1
                    legendHandles = [formula; experiment; band];
                end
            end
        end
        lg = legend(legendHandles, {"formula", experimentLabel, ...
            bandLabel}, "Location", "south", ...
            "Orientation", "horizontal", "Box", "off");
        lg.Layout.Tile = "south";
        lg.FontName = "Times New Roman";
        lg.FontSize = 8;
        switch model
            case "linear", fileStem = "Fig03_Euler_coverage_formula_M30000";
            case "plate", fileStem = "Fig08_plate_coverage_formula_M30000";
            case "nonlinear", fileStem = "Fig10_nonlinear_coverage_formula_M50000";
        end
        exportgraphics(f, fullfile(outputDir, fileStem + ".pdf"), ...
            "ContentType", "vector");
        close(f);
    end
end

manifest = struct();
manifest.schema_version = "nested_M_paper_style_figures_v1";
manifest.source_csv = "data/publication/coverage.csv";
manifest.layout = "two R rows by three confidence-interval method columns";
manifest.band = "pointwise 95% Clopper-Pearson interval centered on empirical coverage";
manifest.display_points = "up to 13 approximately equally spaced response locations";
manifest.generated_at = char(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss"));
manifest.pdf_count = numel(dir(fullfile(outputDir, "*.pdf")));
fid = fopen(fullfile(outputDir, "figure_manifest.json"), "w");
if fid < 0
    error("Cannot write figure manifest.");
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, jsonencode(manifest, PrettyPrint=true), "char");
fprintf("Paper-style nested-M figures written to %s\n", outputDir);
end

function indices = localDisplayIndices(n)
target = min(13, n);
indices = unique(round(linspace(1, n, target))).';
end

function limits = localRowLimits(rows)
values = [rows.coverage; rows.cp95_lower; rows.cp95_upper; ...
    rows.formula_prediction; 0.95];
values = values(isfinite(values));
lo = min(values);
hi = max(values);
padding = max(0.008, 0.06 * (hi - lo));
lo = floor(100 * (lo - padding)) / 100;
hi = ceil(100 * (hi + padding)) / 100;
limits = [lo, hi];
if diff(limits) < 0.04
    center = mean(limits);
    limits = [center - 0.02, center + 0.02];
end
end

function label = localXLabel(model)
if model == "plate"
    label = "response coordinate (m)";
else
    label = "y";
end
end
