function outputs = plot_intervals_pdf_label_only(outputDir, commonPath, platePath)
% Plot recommended-R intervals for all paper examples from archived data.

assert(nargin == 3, 'Use replay_matlab with explicit data and output directories.');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

assert(isfile(commonPath), "Missing Euler/nonlinear data: %s", commonPath);
assert(isfile(platePath), "Missing plate data: %s", platePath);

common = readtable(commonPath, TextType="string", Delimiter=",");
plate = readtable(platePath, TextType="string", Delimiter=",");

specs = [ ...
    localSpec("Euler beam", common(common.model == "linear", :), ...
        "Response y", "Fig05_Euler_recommended_R_intervals_MC_RQMC.pdf"); ...
    localSpec("Kirchhoff plate", plate, "Center deflection w_c (m)", ...
        "Fig09_plate_recommended_R_intervals_MC_RQMC.pdf"); ...
    localSpec("Nonlinear response", common(common.model == "nonlinear", :), ...
        "Response y", "Fig12_nonlinear_recommended_R_intervals_MC_RQMC.pdf") ...
    ];

outputs = strings(numel(specs), 1);
auditRows = table();
for i = 1:numel(specs)
    localValidate(specs(i).data, specs(i).label);
    outputs(i) = localPlot(specs(i), outputDir);
    auditRows = [auditRows; localAuditRows(specs(i))]; %#ok<AGROW>
end
writetable(auditRows, fullfile(outputDir, "response_grid_audit.csv"));

manifest = struct();
manifest.created_at = char(datetime("now", Format="yyyy-MM-dd HH:mm:ss Z"));
manifest.euler_nonlinear_source = "data/publication/intervals_linear_nonlinear.csv";
manifest.plate_source = "data/publication/intervals_plate.csv";
manifest.output_files = arrayfun(@(p) string(localFileName(p)), outputs);
manifest.response_grid_audit = "response_grid_audit.csv";
manifest.rendering = ["MATLAB vector PDF; three method panels; " + ...
    "non-overlapping legends; left and right local enlargements"];
manifest.interpolation = "none";
manifest.validation = "Data interface checked; visual review required";
manifest.viewport = "Density-scale view; interval endpoints outside the view remain unchanged in CSV";
manifest.bootstrap_t_display = "Display only bootstrap-t intervals inside the density viewport and at reference density >=5% of peak; omitted display intervals remain in raw CSV and all coverage trials";
manifestPath = fullfile(outputDir, "plot_manifest.json");
fid = fopen(manifestPath, "w");
assert(fid >= 0, "Cannot write manifest: %s", manifestPath);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, "%s\n", jsonencode(manifest, PrettyPrint=true));
end

function name = localFileName(path)
[~, stem, ext] = fileparts(path);
name = stem + ext;
end

function spec = localSpec(label, data, xLabel, fileName)
spec = struct("label", label, "data", data, "xLabel", xLabel, ...
    "fileName", fileName);
end

function pdfPath = localPlot(spec, outputDir)
methods = ["Student t", "percentile bootstrap", "bootstrap-t"];
methodTitles = ["Student t", "percentile bootstrap", "bootstrap t"];
ink = [34, 39, 44] / 255;
rqmcColor = [31, 92, 151] / 255;
rqmcBand = [164, 198, 225] / 255;
mcColor = [190, 76, 47] / 255;
mcBand = [238, 178, 151] / 255;
gridColor = [224, 228, 231] / 255;

% Preserve the original paper figure; only add two small local views.
fig = figure(Visible="off", Color="white", Units="inches", ...
    Position=[0.4, 0.4, 8.6, 3.1]);
layout = tiledlayout(fig, 1, 3, TileSpacing="compact", Padding="compact");
mainAxes = gobjects(1, 3);
sharedLegendHandles = gobjects(1, 5);
blocks = cell(3, 2);
% Choose the viewing scale from density curves, not extreme CI endpoints.
% This changes only the viewport; all original endpoints remain in the CSV.
densityValues = [double(spec.data.reference_density); double(spec.data.sample_mean)];
commonYLim = [0, 1.18 * max(densityValues(isfinite(densityValues)))];
if spec.label=="Nonlinear response"
    commonYLim(2)=1.55*max(densityValues(isfinite(densityValues)));
end
viewAudit = spec.data(:, {'scheme','method','n','R','y','ci_lower','ci_upper'});
viewAudit.below_view = viewAudit.ci_lower < commonYLim(1);
viewAudit.above_view = viewAudit.ci_upper > commonYLim(2);
writetable(viewAudit, fullfile(outputDir, replace(spec.fileName,'.pdf','_viewport.csv')));

for imethod = 1:numel(methods)
    rqmc = sortrows(spec.data(spec.data.scheme == "rqmc" & ...
        spec.data.method == methods(imethod), :), "y");
    mc = sortrows(spec.data(spec.data.scheme == "mc" & ...
        spec.data.method == methods(imethod), :), "y");
    x = double(rqmc.y).';
    truth = double(rqmc.reference_density).';

    ax = nexttile(layout);
    mainAxes(imethod) = ax;
    blocks{imethod, 1} = rqmc;
    blocks{imethod, 2} = mc;
    hold(ax, "on");
    bandRegion = true(size(truth));
    if methods(imethod)=="bootstrap-t"
        bandRegion = truth >= 0.05*max(truth);
    end
    yUpper = commonYLim(2);

    hTruth = plot(ax, x, truth, Color=ink, LineWidth=1.55, LineStyle="--", ...
        DisplayName="Reference density");
    rqmcRegion=bandRegion;mcRegion=bandRegion;
    if methods(imethod)=="bootstrap-t"
        rqmcRegion=rqmcRegion & double(rqmc.ci_lower).'>=0 & double(rqmc.ci_upper).'<=yUpper;
        mcRegion=mcRegion & double(mc.ci_lower).'>=0 & double(mc.ci_upper).'<=yUpper;
    end
    hRQMCBand = localBand(ax, x, double(rqmc.ci_lower).', ...
        double(rqmc.ci_upper).', rqmcBand, "RQMC 95% interval", ...
        yUpper, rqmcRegion);
    hRQMC = plot(ax, x, double(rqmc.sample_mean).', Color=rqmcColor, ...
        LineWidth=1.15, DisplayName=sprintf("RQMC: n=%d, R=%d", ...
        rqmc.n(1), rqmc.R(1)));
    hMCBand = localBand(ax, x, double(mc.ci_lower).', ...
        double(mc.ci_upper).', mcBand, "MC 95% interval", ...
        yUpper, mcRegion);
    hMC = plot(ax, x, double(mc.sample_mean).', Color=mcColor, ...
        LineWidth=1.15, DisplayName=sprintf("MC: n=%d, R=%d", ...
        mc.n(1), mc.R(1)));
    if imethod == 1
        sharedLegendHandles = [hTruth, hRQMC, hRQMCBand, hMC, hMCBand];
    end

    title(ax, methodTitles(imethod), FontName="Arial", FontSize=10.0, ...
        Color=ink, FontWeight="normal");
    xlabel(ax, spec.xLabel, FontName="Arial", FontSize=9.1, Color=ink);
    if imethod == 1
        ylabel(ax, "pdf", FontName="Arial", FontSize=9.1, Color=ink);
    end
    localStyleAxes(ax, gridColor, ink, 8.4);
    xlim(ax, [min(x), max(x)]);
    ylim(ax, commonYLim);
end

sharedLegend = legend(mainAxes(1), sharedLegendHandles, ...
    {"Reference density", "RQMC mean", "RQMC 95% interval", ...
    "MC mean", "MC 95% interval"}, ...
    Orientation="horizontal", NumColumns=5, Box="off", ...
    FontName="Arial", FontSize=9.1);
sharedLegend.Layout.Tile = "north";
drawnow;

for imethod = 1:numel(methods)
    ax = mainAxes(imethod);
    rqmc = blocks{imethod, 1};
    mc = blocks{imethod, 2};
    x = double(rqmc.y).';
    truth = double(rqmc.reference_density).';
    if methods(imethod)=="bootstrap-t"
        [leftIdx, peakIdx, rightIdx] = localZoomIndices(truth,.08,.15);
    else
        [leftIdx, peakIdx, rightIdx] = localZoomIndices(truth);
    end
    bubbleHeightFraction = 0.43;
    if spec.label == "Nonlinear response"
        bubbleHeightFraction = 0.72;
    end
    localPlotBubble(fig, ax, x, rqmc, mc, leftIdx, "left", ...
        bubbleHeightFraction, ink, rqmcColor, rqmcBand, mcColor, mcBand);
    localPlotBubble(fig, ax, x, rqmc, mc, rightIdx, "right", ...
        bubbleHeightFraction, ink, rqmcColor, rqmcBand, mcColor, mcBand);
    if spec.label == "Nonlinear response"
        localPlotBubble(fig, ax, x, rqmc, mc, peakIdx, "peak", ...
            0.72, ink, rqmcColor, rqmcBand, mcColor, mcBand);
    end
end

pdfPath = string(fullfile(outputDir, spec.fileName));
pngPath = replace(pdfPath, ".pdf", ".png");
exportgraphics(fig, pdfPath, ContentType="vector");
exportgraphics(fig, pngPath, Resolution=300);
close(fig);
end

function localPlotBubble(fig, mainAx, x, rqmc, mc, idx, side, ...
    bubbleHeightFraction, ink, rqmcColor, rqmcBand, mcColor, mcBand)
mainPos = mainAx.Position;
figPixels = getpixelposition(fig);
bubbleWidth = 0.27 * mainPos(3);
if bubbleHeightFraction>=.7,bubbleWidth=.25*mainPos(3);end
bubbleHeight = bubbleWidth * figPixels(3) / figPixels(4);
bubbleY = mainPos(2) + bubbleHeightFraction * mainPos(4);
if side == "left"
    bubbleX = mainPos(1) + 0.025 * mainPos(3);
    sourceIdx = idx(end);
elseif side == "right"
    bubbleX = mainPos(1) + mainPos(3) - bubbleWidth - ...
        0.025 * mainPos(3);
    sourceIdx = idx(1);
else
    bubbleX = mainPos(1) + 0.375 * mainPos(3);
    sourceIdx = idx(ceil(numel(idx) / 2));
end
bubblePos = [bubbleX, bubbleY, bubbleWidth, bubbleHeight];
ax = axes(fig, Position=bubblePos, Color='none');
hold(ax, "on");
rectangle(ax, Position=[0.025, 0.025, 0.95, 0.95], ...
    Curvature=[1, 1], FaceColor="white", EdgeColor="none");
xx = x(idx);
truth = double(rqmc.reference_density(idx)).';
rqmcLower = double(rqmc.ci_lower(idx)).';
rqmcUpper = double(rqmc.ci_upper(idx)).';
mcLower = double(mc.ci_lower(idx)).';
mcUpper = double(mc.ci_upper(idx)).';
if rqmc.method(1)=="bootstrap-t"
    viewport=ylim(mainAx);
    valid=rqmcLower>=viewport(1)&rqmcUpper<=viewport(2);
    rqmcLower(~valid)=NaN;rqmcUpper(~valid)=NaN;
    valid=mcLower>=viewport(1)&mcUpper<=viewport(2);
    mcLower(~valid)=NaN;mcUpper(~valid)=NaN;
end
central = [truth, double(rqmc.sample_mean(idx)).', double(mc.sample_mean(idx)).'];
values = [central,rqmcLower,rqmcUpper,mcLower,mcUpper];
values = values(isfinite(values));
centerBottom = min(central); centerTop = max(central);
centerSpan = max(centerTop-centerBottom, .08*max(abs(central)));
localSpan = max(min(max(values)-min(values),2.2*centerSpan),eps(centerTop));
localBottom = (centerBottom+centerTop-localSpan)/2;
xn = -0.03 + 1.06 * (xx - min(xx)) / (max(xx) - min(xx));
normalizeY = @(v) 0.02 + 0.96 * (v - localBottom) / localSpan;
rqmcLowerN = normalizeY(rqmcLower);
rqmcUpperN = normalizeY(rqmcUpper);
mcLowerN = normalizeY(mcLower);
mcUpperN = normalizeY(mcUpper);

localCircleBand(ax,xn,rqmcLowerN,rqmcUpperN,rqmcBand);
localCircleBand(ax,xn,mcLowerN,mcUpperN,mcBand);
localCircleLine(ax,xn,normalizeY(truth),ink,.85,'--');
localCircleLine(ax,xn,normalizeY(double(rqmc.sample_mean(idx)).'),rqmcColor,.72,'-');
localCircleLine(ax,xn,normalizeY(double(mc.sample_mean(idx)).'),mcColor,.72,'-');
rectangle(ax, Position=[0.025, 0.025, 0.95, 0.95], ...
    Curvature=[1, 1], FaceColor="none", ...
    EdgeColor=[0.35, 0.38, 0.40], LineWidth=0.75);
axis(ax, [0, 1, 0, 1]);
axis(ax, "off");

mainXLim = xlim(mainAx);
mainYLim = ylim(mainAx);
sourceX = mainPos(1) + ...
    (x(sourceIdx) - mainXLim(1)) / diff(mainXLim) * mainPos(3);
sourceValue = double(rqmc.reference_density(sourceIdx));
sourceY = mainPos(2) + ...
    (sourceValue - mainYLim(1)) / diff(mainYLim) * mainPos(4);
center = [bubbleX + bubbleWidth / 2, bubbleY + bubbleHeight / 2];
delta = [sourceX - center(1), sourceY - center(2)];
deltaPixels = [delta(1) * figPixels(3), delta(2) * figPixels(4)];
distancePixels = hypot(deltaPixels(1), deltaPixels(2));
if distancePixels > 0
    radiusPixels = 0.475 * bubbleWidth * figPixels(3);
    boundaryPixels = radiusPixels * deltaPixels / distancePixels;
    start = center + [boundaryPixels(1) / figPixels(3), ...
        boundaryPixels(2) / figPixels(4)];
else
    start = center;
end
annotation(fig, "arrow", [start(1), sourceX], [start(2), sourceY], ...
    Color=[0.35, 0.38, 0.40], LineWidth=0.65, HeadLength=5, ...
    HeadWidth=5);
end

function localCircleBand(ax,x,lower,upper,color)
angle=(0:179)*2*pi/180;
disc=polyshape(.5+.475*cos(angle),.5+.475*sin(angle));
valid=isfinite(lower)&isfinite(upper);
edges=diff([false,valid,false]);starts=find(edges==1);stops=find(edges==-1)-1;
for k=1:numel(starts)
    ii=starts(k):stops(k);if numel(ii)<2,continue;end
    band=polyshape([x(ii),fliplr(x(ii))],[lower(ii),fliplr(upper(ii))]);
    clipped=intersect(band,disc);
    plot(ax,clipped,FaceColor=color,FaceAlpha=.48,EdgeColor='none',HandleVisibility='off');
end
end

function localCircleLine(ax,x,y,color,width,style)
xx=linspace(min(x),max(x),801);yy=interp1(x,y,xx,'linear');
yy((xx-.5).^2+(yy-.5).^2>.475^2)=NaN;
plot(ax,xx,yy,Color=color,LineWidth=width,LineStyle=style);
end

function [leftIdx, peakZoomIdx, rightIdx] = localZoomIndices(truth,lowerFraction,upperFraction)
if nargin<2,lowerFraction=.015;upperFraction=.045;end
truth = truth(:).';
[peak, peakIdx] = max(truth);
leftCandidates = find((1:numel(truth)) < peakIdx & ...
    truth >= lowerFraction * peak & truth <= upperFraction * peak);
rightCandidates = find((1:numel(truth)) > peakIdx & ...
    truth >= lowerFraction * peak & truth <= upperFraction * peak);
if numel(leftCandidates) < 5
    leftCandidates = find((1:numel(truth)) < peakIdx & ...
        truth >= lowerFraction*peak & truth <= max(.25,upperFraction)*peak);
end
if numel(rightCandidates) < 5
    rightCandidates = find((1:numel(truth)) > peakIdx & ...
        truth >= lowerFraction*peak & truth <= max(.25,upperFraction)*peak);
end
assert(numel(leftCandidates)>=2&&numel(rightCandidates)>=2,'Response grid does not resolve both enlargement regions');
leftIdx = min(leftCandidates):max(leftCandidates);
rightIdx = min(rightCandidates):max(rightCandidates);
peakHalfWidth = max(4, round(0.025 * numel(truth)));
peakZoomIdx = max(1, peakIdx - peakHalfWidth): ...
    min(numel(truth), peakIdx + peakHalfWidth);
end

function h = localBand(ax, x, lower, upper, color, displayName, yUpper, region)
valid = isfinite(lower) & isfinite(upper) & region;
edges = diff([false, valid, false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
if isempty(starts)
    h=fill(ax,NaN,NaN,color,EdgeColor='none',FaceAlpha=.30,DisplayName=displayName);
    return
end
h = gobjects(1);
for i = 1:numel(starts)
    idx = starts(i):stops(i);
    args = {"EdgeColor", "none", "FaceAlpha", 0.30};
    if i == 1 && strlength(string(displayName)) > 0
        args = [args, {"DisplayName", displayName}]; %#ok<AGROW>
    else
        args = [args, {"HandleVisibility", "off"}]; %#ok<AGROW>
    end
    patch = fill(ax, [x(idx), fliplr(x(idx))], ...
        [lower(idx), fliplr(upper(idx))], color, args{:});
    if i == 1
        h = patch;
    end
end
end

function localStyleAxes(ax, gridColor, ink, fontSize)
grid(ax, "on");
ax.GridColor = gridColor;
ax.GridAlpha = 0.52;
ax.MinorGridAlpha = 0.16;
ax.FontName = "Arial";
ax.FontSize = fontSize;
ax.XColor = ink;
ax.YColor = ink;
ax.Box = "on";
end

function localValidate(T, label)
required = ["scheme", "method", "n", "R", "B", "y", ...
    "reference_density", "sample_mean", "ci_lower", "ci_upper"];
assert(all(ismember(required, string(T.Properties.VariableNames))), ...
    "%s data are missing required columns.", label);
methods = ["Student t", "percentile bootstrap", "bootstrap-t"];
for scheme = ["rqmc", "mc"]
    for method = methods
        rows = T(T.scheme == scheme & T.method == method, :);
        assert(height(rows) == 201, ...
            "%s %s %s must contain 201 response points.", ...
            label, scheme, method);
        y = sort(double(rows.y));
        assert(all(diff(y) > 0), "%s response grid is not unique.", label);
    end
end
end

function rows = localAuditRows(spec)
methods = ["Student t", "percentile bootstrap", "bootstrap-t"];
rows = table();
for scheme = ["rqmc", "mc"]
    for method = methods
        block = spec.data(spec.data.scheme == scheme & ...
            spec.data.method == method, :);
        row = table(string(spec.label), string(scheme), string(method), ...
            double(height(block)), double(numel(unique(block.y))), ...
            double(min(block.y)), double(max(block.y)), ...
            VariableNames=["example", "scheme", "method", ...
            "row_count", "unique_response_points", "response_min", ...
            "response_max"]);
        rows = [rows; row]; %#ok<AGROW>
    end
end
end
