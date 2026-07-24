package cra.autosplit;

import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.plugin.filter.ThresholdToSelection;
import ij.process.BinaryProcessor;
import ij.process.ByteProcessor;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/** Marker-controlled geodesic partitioning for CRA lobule candidates. */
public final class AutoSplitEngine {
    public static final String AUTO_PREFIX = "AUTO_REVIEW_";

    public static final class Seed {
        public final String region;
        public final double x;
        public final double y;

        public Seed(String region, double x, double y) {
            this.region = region;
            this.x = x;
            this.y = y;
        }
    }

    public static final class Options {
        public int barrierWidthPixels = 12;
        public int minimumComponentPixels = 80;
        public boolean createPL = true;
        public boolean createML = true;
        public boolean createGL = true;
        public boolean createDerivedTotals = true;
    }

    public static final class Result {
        public final List<Roi> candidates = new ArrayList<Roi>();
        public final List<String> warnings = new ArrayList<String>();
        public final List<String> diagnostics = new ArrayList<String>();
        public Rectangle crop;
        public int gridWidth;
        public int gridHeight;
        public long elapsedMillis;
        public long cbPixels;
        public long assignedCbPixels;
        public long plMlOverlapPixels;
        public long plGlOverlapPixels;
        public long mlGlOverlapPixels;
        public long glwmOutsideCbPixels;
        public long wmOutsideGlwmPixels;
        public long plOutsideCbPixels;
        public long plWmOverlapPixels;

        public double assignedPercent() {
            if (cbPixels == 0) return 0.0;
            return 100.0 * assignedCbPixels / cbPixels;
        }
    }

    /* Preview ROIs keep their AUTO_REVIEW_ names, but use the established CRA
     * layer colours.  A single purple colour made adjacent PL/ML/GL outlines
     * look as if the masks overlapped even when the raster masks were disjoint. */
    private static final Color PREVIEW_PL_COLOR = new Color(255, 64, 96);
    private static final Color PREVIEW_ML_COLOR = new Color(0, 210, 90);
    private static final Color PREVIEW_GL_COLOR = new Color(255, 190, 0);
    private static final Color PREVIEW_TOTAL_COLOR = new Color(0, 190, 255);

    private AutoSplitEngine() {}

    public static Result generate(Map<String, Roi> namedRois, List<Seed> seeds, Options options) {
        long started = System.currentTimeMillis();
        requireArea(namedRois, "Cb");
        requireArea(namedRois, "WM");
        requireArea(namedRois, "GL+WM");
        Roi plTotal = require(namedRois, "PL_TOTAL");
        if (seeds == null || seeds.size() < 2)
            throw new IllegalArgumentException("At least two valid single-point lobule seeds are required.");
        validateSeeds(seeds);

        Options safe = sanitize(options);
        Roi cb = namedRois.get("Cb");
        Roi wm = namedRois.get("WM");
        Roi glwm = namedRois.get("GL+WM");
        Rectangle crop = paddedBounds(cb.getBounds(), 2);
        int width = crop.width;
        int height = crop.height;

        byte[] cbMask = roiToMask(cb, crop, 1);
        byte[] rawWmMask = roiToMask(wm, crop, 1);
        byte[] rawGlwmMask = roiToMask(glwm, crop, 1);
        byte[] rawPlMask = roiToMask(plTotal, crop, plTotal.isArea() ? 1 : 3);
        Roi suppliedMlTotal = namedRois.get("ML_TOTAL");
        byte[] rawSuppliedMlMask = suppliedMlTotal == null ? null :
            roiToMask(suppliedMlTotal, crop, suppliedMlTotal.isArea() ? 1 : 3);

        Result result = new Result();
        result.crop = crop;
        result.gridWidth = width;
        result.gridHeight = height;
        result.glwmOutsideCbPixels = countMask(subtractMask(rawGlwmMask, cbMask));
        result.wmOutsideGlwmPixels = countMask(subtractMask(rawWmMask, rawGlwmMask));
        result.plOutsideCbPixels = countMask(subtractMask(rawPlMask, cbMask));
        result.plWmOverlapPixels = countOverlap(rawPlMask, rawWmMask);
        auditInputNesting("GL+WM outside Cb", result.glwmOutsideCbPixels,
            countMask(rawGlwmMask), result);
        auditInputNesting("WM outside GL+WM", result.wmOutsideGlwmPixels,
            countMask(rawWmMask), result);
        auditInputNesting("PL_TOTAL outside Cb", result.plOutsideCbPixels,
            countMask(rawPlMask), result);
        auditInputNesting("PL_TOTAL overlapping WM", result.plWmOverlapPixels,
            countMask(rawPlMask), result);

        /* Enforce the anatomical nesting before subtraction.  This is intentionally
         * not XOR: small drawing overshoots must never become false layer islands. */
        byte[] glwmMask = andMask(cbMask, rawGlwmMask);
        byte[] wmMask = andMask(glwmMask, rawWmMask);
        byte[] plMask = andMask(cbMask, rawPlMask);

        /* PL is a measured layer of its own, not part of ML or GL.  ML_TOTAL is
         * interpreted as the molecular-side area between the pial/Cb boundary
         * and PL_TOTAL: the pial-side band only, never the thin "inner ML" sliver
         * trapped between PL and the GL+WM boundary.  PL_TOTAL is traced as ONE
         * continuous open curve across the whole cortical ribbon (it runs through
         * every fissure rather than closing into a separate loop per lobule), so
         * naive hole-filling can never find a lobule's interior pocket: an open
         * arc does not separate the plane, so a border-seeded flood can always
         * sneak around through the trace's two far-away open endpoints and reach
         * every pocket.  Seal those two endpoints to GL+WM instead, then keep
         * only the candidate ML pixels reachable from the pial (outer Cb) surface
         * without crossing PL.  This works regardless of whether ML_TOTAL is
         * supplied or derived, and regardless of how many lobules PL winds through. */
        byte[] mlCandidate;
        if (rawSuppliedMlMask != null) {
            byte[] suppliedInsideCb = andMask(cbMask, rawSuppliedMlMask);
            mlCandidate = subtractMask(subtractMask(suppliedInsideCb, plMask), glwmMask);
            result.diagnostics.add("ML mask source: supplied ML_TOTAL (clipped to Cb and made exclusive from PL/GL+WM)");
        } else {
            mlCandidate = subtractMask(subtractMask(cbMask, glwmMask), plMask);
            result.diagnostics.add("ML mask source: PL_TOTAL exterior envelope (pial-reachable tissue only)");
        }
        byte[] mlMask = pialSideMlMask(cbMask, mlCandidate, glwmMask, plTotal, crop);
        long innerMlSliverPixels = countMask(mlCandidate) - countMask(mlMask);
        if (innerMlSliverPixels > 0)
            result.diagnostics.add("Inner-ML sliver removed: " + innerMlSliverPixels +
                " px (kept pial-side ML only).");
        byte[] mlTopologyMask = mlMask;
        byte[] glTopologyMask = subtractMask(glwmMask, wmMask);
        byte[] glMask = subtractMask(glTopologyMask, plMask);
        long plMlOverlap = countOverlap(plMask, mlMask);
        long plGlOverlap = countOverlap(plMask, glMask);
        long mlGlOverlap = countOverlap(mlMask, glMask);
        if (plMlOverlap != 0 || plGlOverlap != 0 || mlGlOverlap != 0)
            throw new IllegalStateException("Internal layer-mask overlap detected; no candidates were committed.");
        byte[] barriers = separatorMask(namedRois, crop, safe.barrierWidthPixels);
        int separatorCount = countSeparators(namedRois);

        /* The parent components keep the detached inferior tissue from borrowing
         * the name of a geometrically close lobule in the main tissue. */
        Labeling parentLabels = labelComponents(glwmMask, null, width, height, 500, seeds, crop);
        int[] seedParents = assignSeedParents(seeds, parentLabels, crop, width, height);

        if (separatorCount == 0)
            result.warnings.add("No FISSURE_ or DIVIDER_ guide was found; boundaries use seed distance only.");
        long layerPixels = 0;
        if (safe.createPL) layerPixels += countMask(plMask);
        if (safe.createML) layerPixels += countMask(mlMask);
        if (safe.createGL) layerPixels += countMask(glMask);
        result.cbPixels = layerPixels;
        result.plMlOverlapPixels = countOverlap(plMask, mlMask);
        result.plGlOverlapPixels = countOverlap(plMask, glMask);
        result.mlGlOverlapPixels = countOverlap(mlMask, glMask);
        result.diagnostics.add("Layer-mask overlap audit: PL-ML=" + result.plMlOverlapPixels +
            ", PL-GL=" + result.plGlOverlapPixels +
            ", ML-GL=" + result.mlGlOverlapPixels + " pixels");
        result.diagnostics.add("Layer overlap check: PL∩ML=0, PL∩GL=0, ML∩GL=0 pixels");

        if (safe.createDerivedTotals) {
            if (!namedRois.containsKey("ML_TOTAL")) {
                addDerivedTotal("ML_TOTAL", width, height, crop, new PixelRule() {
                    public boolean include(int p) {
                        return (mlMask[p] & 0xff) != 0;
                    }
                }, result.candidates);
            }
            if (!namedRois.containsKey("GL_TOTAL")) {
                addDerivedTotal("GL_TOTAL", width, height, crop, new PixelRule() {
                    public boolean include(int p) {
                        return (glMask[p] & 0xff) != 0;
                    }
                }, result.candidates);
            }
        }

        long assigned = 0;
        if (safe.createPL) {
            if (plTotal.isLine())
                assigned += splitPLDirect(plTotal, plMask, barriers, seeds, parentLabels,
                    seedParents, crop, Math.max(8, safe.minimumComponentPixels / 8),
                    safe.barrierWidthPixels, result);
            else
                assigned += splitLayer("PL", plMask, plMask, barriers, seeds, parentLabels,
                    seedParents, crop, Math.max(8, safe.minimumComponentPixels / 8),
                    safe.barrierWidthPixels, result);
        }
        if (safe.createML)
            assigned += splitLayer("ML", mlTopologyMask, mlMask, barriers, seeds, parentLabels,
                seedParents, crop, safe.minimumComponentPixels,
                safe.barrierWidthPixels, result);
        if (safe.createGL)
            assigned += splitLayer("GL", glTopologyMask, glMask, barriers, seeds, parentLabels,
                seedParents, crop, safe.minimumComponentPixels,
                safe.barrierWidthPixels, result);
        result.assignedCbPixels = assigned;
        if (result.assignedPercent() < 99.0)
            result.warnings.add(String.format("Meaningful layer components cover %.2f%% of PL+ML+GL pixels; tiny artifacts/barrier pixels were excluded.", result.assignedPercent()));

        result.elapsedMillis = System.currentTimeMillis() - started;
        if (result.candidates.isEmpty())
            result.warnings.add("No candidate ROI was generated.");
        return result;
    }

    private static void auditInputNesting(String label, long violation, long denominator,
            Result result) {
        double percent = denominator <= 0 ? 0.0 : 100.0 * violation / denominator;
        result.diagnostics.add(label + ": " + violation + " px (" +
            String.format("%.3f%%", percent) + ")");
        if (violation > 1000 && percent > 10.0)
            throw new IllegalArgumentException(label + " is too large: " + violation +
                " pixels (" + String.format("%.2f%%", percent) +
                "). Fix the whole/total ROI nesting before Preview.");
        if (violation > 100 && percent > 0.5)
            result.warnings.add(label + " was clipped: " + violation + " px (" +
                String.format("%.2f%%", percent) + "). Inspect the corresponding total ROI.");
    }

    private static Options sanitize(Options source) {
        Options out = new Options();
        if (source != null) {
            out.barrierWidthPixels = Math.max(1, Math.min(200, source.barrierWidthPixels));
            out.minimumComponentPixels = Math.max(1, source.minimumComponentPixels);
            out.createPL = source.createPL;
            out.createML = source.createML;
            out.createGL = source.createGL;
            out.createDerivedTotals = source.createDerivedTotals;
        }
        return out;
    }

    private static void validateSeeds(List<Seed> seeds) {
        for (int i = 0; i < seeds.size(); i++) {
            Seed a = seeds.get(i);
            if (a.region == null || a.region.trim().isEmpty())
                throw new IllegalArgumentException("A seed has no lobule name.");
            for (int j = i + 1; j < seeds.size(); j++) {
                Seed b = seeds.get(j);
                if (a.region.equals(b.region))
                    throw new IllegalArgumentException("Duplicate seed region: " + a.region);
                double ddx = a.x - b.x;
                double ddy = a.y - b.y;
                if (ddx * ddx + ddy * ddy < 9.0)
                    throw new IllegalArgumentException("Seeds " + a.region + " and " + b.region + " overlap.");
            }
        }
    }

    private static Roi require(Map<String, Roi> rois, String name) {
        Roi roi = rois.get(name);
        if (roi == null) throw new IllegalArgumentException("Required ROI is missing: " + name);
        return roi;
    }

    private static Roi requireArea(Map<String, Roi> rois, String name) {
        Roi roi = require(rois, name);
        if (!roi.isArea()) throw new IllegalArgumentException(name + " must be an area ROI.");
        return roi;
    }

    private static Rectangle paddedBounds(Rectangle source, int padding) {
        return new Rectangle(Math.max(0, source.x - padding), Math.max(0, source.y - padding),
            source.width + padding * 2, source.height + padding * 2);
    }

    private static byte[] roiToMask(Roi source, Rectangle crop, int lineWidth) {
        ByteProcessor processor = new ByteProcessor(crop.width, crop.height);
        processor.setValue(255);
        Roi local = (Roi) source.clone();
        Rectangle bounds = local.getBounds();
        local.setLocation(bounds.x - crop.x, bounds.y - crop.y);
        if (source.isArea()) {
            processor.fill(local);
        } else {
            processor.setColor(255);
            processor.setLineWidth(Math.max(1, lineWidth));
            processor.draw(local);
        }
        return (byte[]) processor.getPixels();
    }

    private static byte[] andMask(byte[] left, byte[] right) {
        byte[] out = new byte[left.length];
        for (int i = 0; i < out.length; i++) {
            if ((left[i] & 0xff) != 0 && (right[i] & 0xff) != 0) out[i] = (byte) 255;
        }
        return out;
    }

    private static byte[] subtractMask(byte[] left, byte[] right) {
        byte[] out = new byte[left.length];
        for (int i = 0; i < out.length; i++) {
            if ((left[i] & 0xff) != 0 && (right[i] & 0xff) == 0) out[i] = (byte) 255;
        }
        return out;
    }

    /** Keep only the candidate ML pixels that can reach the pial (outer Cb)
     * surface without crossing PL.  PL_TOTAL traces continuously through every
     * fissure as a single open curve, so its two true endpoints are sealed to
     * the nearest GL+WM pixel first (an open arc alone can never wall off a
     * lobule's interior); the continuous PL stroke, already removed from
     * mlCandidate, blocks the rest of the boundary.  A 4-connected flood from
     * the pial surface then reaches the outer molecular band and leaves any
     * inner-ML sliver trapped against GL+WM behind. */
    private static byte[] pialSideMlMask(byte[] cbMask, byte[] mlCandidate, byte[] glwmMask,
            Roi plTotal, Rectangle crop) {
        int width = crop.width;
        int height = crop.height;
        byte[] blocker = new byte[mlCandidate.length];
        FloatPolygon points = plTotal.getFloatPolygon();
        if (points != null && points.npoints >= 2) {
            sealEndpointToGlwm(blocker, glwmMask, width, height,
                points.xpoints[0] - crop.x, points.ypoints[0] - crop.y);
            sealEndpointToGlwm(blocker, glwmMask, width, height,
                points.xpoints[points.npoints - 1] - crop.x,
                points.ypoints[points.npoints - 1] - crop.y);
        }
        boolean[] reached = new boolean[mlCandidate.length];
        int[] queue = new int[mlCandidate.length];
        int head = 0;
        int tail = 0;
        for (int y = 0; y < height; y++) {
            int row = y * width;
            for (int x = 0; x < width; x++) {
                int p = row + x;
                if ((mlCandidate[p] & 0xff) == 0 || (blocker[p] & 0xff) != 0 || reached[p]) continue;
                if (!touchesBackground(cbMask, x, y, width, height)) continue;
                reached[p] = true;
                queue[tail++] = p;
            }
        }
        final int[] nx = {-1, 1, 0, 0};
        final int[] ny = {0, 0, -1, 1};
        while (head < tail) {
            int p = queue[head++];
            int x = p % width;
            int y = p / width;
            for (int d = 0; d < 4; d++) {
                int xx = x + nx[d];
                int yy = y + ny[d];
                if (xx < 0 || yy < 0 || xx >= width || yy >= height) continue;
                int q = yy * width + xx;
                if (reached[q] || (mlCandidate[q] & 0xff) == 0 || (blocker[q] & 0xff) != 0) continue;
                reached[q] = true;
                queue[tail++] = q;
            }
        }
        byte[] out = new byte[mlCandidate.length];
        for (int i = 0; i < out.length; i++) if (reached[i]) out[i] = (byte) 255;
        return out;
    }

    private static boolean touchesBackground(byte[] cbMask, int x, int y, int width, int height) {
        if (x == 0 || y == 0 || x == width - 1 || y == height - 1) return true;
        return (cbMask[y * width + (x - 1)] & 0xff) == 0
            || (cbMask[y * width + (x + 1)] & 0xff) == 0
            || (cbMask[(y - 1) * width + x] & 0xff) == 0
            || (cbMask[(y + 1) * width + x] & 0xff) == 0;
    }

    /** Draw a short barrier from an open PL endpoint to the nearest GL+WM pixel so
     * the inner-ML sliver is closed off before the pial flood fill. */
    private static void sealEndpointToGlwm(byte[] blocker, byte[] glwmMask,
            int width, int height, double endpointX, double endpointY) {
        int sx = (int) Math.round(endpointX);
        int sy = (int) Math.round(endpointY);
        sx = Math.max(0, Math.min(width - 1, sx));
        sy = Math.max(0, Math.min(height - 1, sy));
        int target = -1;
        int limit = Math.max(width, height);
        for (int radius = 0; radius <= limit && target < 0; radius++) {
            for (int yy = sy - radius; yy <= sy + radius && target < 0; yy++) {
                if (yy < 0 || yy >= height) continue;
                for (int xx = sx - radius; xx <= sx + radius; xx++) {
                    if (xx < 0 || xx >= width) continue;
                    if (radius > 0 && Math.abs(xx - sx) != radius && Math.abs(yy - sy) != radius) continue;
                    if ((glwmMask[yy * width + xx] & 0xff) != 0) {
                        target = yy * width + xx;
                        break;
                    }
                }
            }
        }
        if (target < 0) return;
        drawSealSegment(blocker, width, height, sx, sy, target % width, target / width, 2);
    }

    private static void drawSealSegment(byte[] mask, int width, int height,
            int x0, int y0, int x1, int y1, int radius) {
        int steps = Math.max(Math.abs(x1 - x0), Math.abs(y1 - y0));
        if (steps == 0) steps = 1;
        for (int step = 0; step <= steps; step++) {
            double t = step / (double) steps;
            int cx = (int) Math.round(x0 + t * (x1 - x0));
            int cy = (int) Math.round(y0 + t * (y1 - y0));
            for (int oy = -radius; oy <= radius; oy++) {
                for (int ox = -radius; ox <= radius; ox++) {
                    if (ox * ox + oy * oy > radius * radius) continue;
                    int xx = cx + ox;
                    int yy = cy + oy;
                    if (xx >= 0 && yy >= 0 && xx < width && yy < height)
                        mask[yy * width + xx] = (byte) 255;
                }
            }
        }
    }

    private static long countMask(byte[] mask) {
        long count = 0;
        for (byte value : mask) if ((value & 0xff) != 0) count++;
        return count;
    }

    private static long countOverlap(byte[] left, byte[] right) {
        long count = 0;
        for (int i = 0; i < left.length; i++) {
            if ((left[i] & 0xff) != 0 && (right[i] & 0xff) != 0) count++;
        }
        return count;
    }

    private static void applyPreviewStyle(Roi roi, String type) {
        roi.setFillColor(null);
        if ("PL".equals(type)) {
            roi.setStrokeColor(PREVIEW_PL_COLOR);
            roi.setStrokeWidth(3.0);
        } else if ("ML".equals(type)) {
            roi.setStrokeColor(PREVIEW_ML_COLOR);
            roi.setStrokeWidth(1.25);
        } else if ("GL".equals(type)) {
            roi.setStrokeColor(PREVIEW_GL_COLOR);
            roi.setStrokeWidth(1.25);
        } else {
            roi.setStrokeColor(PREVIEW_TOTAL_COLOR);
            roi.setStrokeWidth(1.5);
        }
    }

    private static int countSeparators(Map<String, Roi> namedRois) {
        int count = 0;
        for (Map.Entry<String, Roi> entry : namedRois.entrySet()) {
            String name = entry.getKey();
            if ((name.startsWith("FISSURE_") || name.startsWith("DIVIDER_")) && entry.getValue().isLine()) count++;
        }
        return count;
    }

    private static byte[] separatorMask(Map<String, Roi> namedRois, Rectangle crop, int widthPixels) {
        byte[] out = new byte[crop.width * crop.height];
        for (Map.Entry<String, Roi> entry : namedRois.entrySet()) {
            String name = entry.getKey();
            if (!name.startsWith("FISSURE_") && !name.startsWith("DIVIDER_")) continue;
            Roi roi = entry.getValue();
            if (!roi.isLine()) continue;
            byte[] one = roiToMask(roi, crop, widthPixels);
            for (int i = 0; i < out.length; i++) if ((one[i] & 0xff) != 0) out[i] = (byte) 255;
        }
        return out;
    }

    private static final class Component {
        int id;
        int area;
        int minX = Integer.MAX_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxX = -1;
        int maxY = -1;
        long sumX;
        long sumY;
        double[] minSeedDistanceSquared;

        double centerX() { return area == 0 ? 0.0 : sumX / (double) area; }
        double centerY() { return area == 0 ? 0.0 : sumY / (double) area; }
    }

    private static final class Labeling {
        final int[] labels;
        final List<Component> components;

        Labeling(int[] labels, List<Component> components) {
            this.labels = labels;
            this.components = components;
        }
    }

    /** Full-resolution connected components.  Barriers are removed only for the
     * layer currently being split, so an ML-only guide cannot corrupt GL or PL. */
    private static Labeling labelComponents(byte[] mask, byte[] barriers, int width, int height,
            int minimumArea, List<Seed> seeds, Rectangle crop) {
        int[] labels = new int[mask.length];
        int[] queue = new int[mask.length];
        List<Component> all = new ArrayList<Component>();
        int nextId = 0;
        final int[] nx = {-1, 0, 1, -1, 1, -1, 0, 1};
        final int[] ny = {-1, -1, -1, 0, 0, 1, 1, 1};

        for (int start = 0; start < mask.length; start++) {
            if ((mask[start] & 0xff) == 0 || labels[start] != 0 ||
                    (barriers != null && (barriers[start] & 0xff) != 0)) continue;
            Component component = new Component();
            component.id = ++nextId;
            if (seeds != null) {
                component.minSeedDistanceSquared = new double[seeds.size()];
                Arrays.fill(component.minSeedDistanceSquared, Double.POSITIVE_INFINITY);
            }
            int head = 0;
            int tail = 0;
            queue[tail++] = start;
            labels[start] = component.id;
            while (head < tail) {
                int pixel = queue[head++];
                int x = pixel % width;
                int y = pixel / width;
                component.area++;
                component.sumX += x;
                component.sumY += y;
                if (x < component.minX) component.minX = x;
                if (x > component.maxX) component.maxX = x;
                if (y < component.minY) component.minY = y;
                if (y > component.maxY) component.maxY = y;
                if (seeds != null) {
                    double imageX = x + crop.x;
                    double imageY = y + crop.y;
                    for (int s = 0; s < seeds.size(); s++) {
                        double dx = imageX - seeds.get(s).x;
                        double dy = imageY - seeds.get(s).y;
                        double distance = dx * dx + dy * dy;
                        if (distance < component.minSeedDistanceSquared[s])
                            component.minSeedDistanceSquared[s] = distance;
                    }
                }
                for (int d = 0; d < 8; d++) {
                    int xx = x + nx[d];
                    int yy = y + ny[d];
                    if (xx < 0 || yy < 0 || xx >= width || yy >= height) continue;
                    int adjacent = yy * width + xx;
                    if (labels[adjacent] != 0 || (mask[adjacent] & 0xff) == 0 ||
                            (barriers != null && (barriers[adjacent] & 0xff) != 0)) continue;
                    labels[adjacent] = component.id;
                    queue[tail++] = adjacent;
                }
            }
            if (component.area >= minimumArea) all.add(component);
        }
        return new Labeling(labels, all);
    }

    private static int[] assignSeedParents(List<Seed> seeds, Labeling parentLabels, Rectangle crop,
            int width, int height) {
        int[] assignment = new int[seeds.size()];
        Arrays.fill(assignment, -1);
        List<Component> parents = parentLabels.components;
        if (parents.isEmpty()) return assignment;
        for (int s = 0; s < seeds.size(); s++) {
            double best = Double.POSITIVE_INFINITY;
            for (int p = 0; p < parents.size(); p++) {
                Component parent = parents.get(p);
                double distance = parent.minSeedDistanceSquared == null
                    ? Double.POSITIVE_INFINITY : parent.minSeedDistanceSquared[s];
                if (distance < best) {
                    best = distance;
                    assignment[s] = p;
                }
            }
        }
        return assignment;
    }

    private static int componentParent(Component component, Labeling layer, Labeling parentLabels,
            int width, int height) {
        List<Component> parents = parentLabels.components;
        if (parents.isEmpty()) return -1;
        int[] idToIndex = new int[parents.size() + 1];
        Arrays.fill(idToIndex, -1);
        int maximumId = 0;
        for (Component parent : parents) if (parent.id > maximumId) maximumId = parent.id;
        if (maximumId >= idToIndex.length) {
            idToIndex = new int[maximumId + 1];
            Arrays.fill(idToIndex, -1);
        }
        for (int i = 0; i < parents.size(); i++) idToIndex[parents.get(i).id] = i;
        int[] contacts = new int[parents.size()];
        final int radius = 24;
        final int stride = component.area > 1000000 ? 5 : (component.area > 200000 ? 3 : 1);
        for (int y = component.minY; y <= component.maxY; y += stride) {
            int row = y * width;
            for (int x = component.minX; x <= component.maxX; x += stride) {
                if (layer.labels[row + x] != component.id) continue;
                for (int distance = 0; distance <= radius; distance += 4) {
                    incrementParentContact(contacts, idToIndex, parentLabels.labels,
                        x - distance, y, width, height);
                    incrementParentContact(contacts, idToIndex, parentLabels.labels,
                        x + distance, y, width, height);
                    incrementParentContact(contacts, idToIndex, parentLabels.labels,
                        x, y - distance, width, height);
                    incrementParentContact(contacts, idToIndex, parentLabels.labels,
                        x, y + distance, width, height);
                }
            }
        }
        int contactWinner = -1;
        int mostContacts = 0;
        for (int i = 0; i < contacts.length; i++) {
            if (contacts[i] > mostContacts) {
                mostContacts = contacts[i];
                contactWinner = i;
            }
        }
        if (contactWinner >= 0) return contactWinner;

        int bestIndex = -1;
        double best = Double.POSITIVE_INFINITY;
        for (int p = 0; p < parents.size(); p++) {
            Component parent = parents.get(p);
            double dx = component.centerX() - parent.centerX();
            double dy = component.centerY() - parent.centerY();
            double distance = dx * dx + dy * dy;
            if (distance < best) {
                best = distance;
                bestIndex = p;
            }
        }
        return bestIndex;
    }

    private static void incrementParentContact(int[] contacts, int[] idToIndex,
            int[] parentPixels, int x, int y, int width, int height) {
        if (x < 0 || y < 0 || x >= width || y >= height) return;
        int parentId = parentPixels[y * width + x];
        if (parentId <= 0 || parentId >= idToIndex.length) return;
        int index = idToIndex[parentId];
        if (index >= 0) contacts[index]++;
    }

    private static long splitLayer(String type, byte[] topologyMask, byte[] outputMask, byte[] barriers,
            List<Seed> seeds, Labeling parentLabels, int[] seedParents, Rectangle crop,
            int requestedMinimumArea, int barrierWidthPixels, Result result) {
        int width = crop.width;
        int height = crop.height;
        long layerArea = countMask(topologyMask);
        int dynamicMinimum = (int) Math.max(requestedMinimumArea, Math.round(layerArea * 0.0001));
        Labeling labeling = labelComponents(topologyMask, barriers, width, height, dynamicMinimum, seeds, crop);
        int maximumTopologyId = 0;
        for (Component component : labeling.components)
            if (component.id > maximumTopologyId) maximumTopologyId = component.id;
        int[] ownerByTopologyId = new int[maximumTopologyId + 1];
        Arrays.fill(ownerByTopologyId, -1);

        for (int c = 0; c < labeling.components.size(); c++) {
            Component component = labeling.components.get(c);
            int parent = componentParent(component, labeling, parentLabels, width, height);
            int bestSeed = -1;
            double best = Double.POSITIVE_INFINITY;
            double second = Double.POSITIVE_INFINITY;
            for (int s = 0; s < seeds.size(); s++) {
                if (parent >= 0 && seedParents[s] >= 0 && parent != seedParents[s]) continue;
                double distance = component.minSeedDistanceSquared[s];
                if (distance < best) {
                    second = best;
                    best = distance;
                    bestSeed = s;
                } else if (distance < second) {
                    second = distance;
                }
            }
            if (bestSeed < 0) {
                result.warnings.add(type + " component " + component.id +
                    " could not be assigned within its tissue parent.");
                continue;
            }
            if (component.id < ownerByTopologyId.length)
                ownerByTopologyId[component.id] = bestSeed;
            result.diagnostics.add(type + " component " + component.id + ": area=" + component.area +
                ", parent=" + parent + ", owner=" + seeds.get(bestSeed).region +
                ", nearestDistance=" + String.format("%.1f px", Math.sqrt(best)));
            if (second < Double.POSITIVE_INFINITY) {
                double margin = Math.sqrt(second) - Math.sqrt(best);
                if (margin < 20.0)
                    result.warnings.add(type + " component " + component.id +
                        " has a low seed-distance margin (" + String.format("%.1f px", margin) +
                        "); inspect its boundary carefully.");
            }
        }

        /* Guides are assignment barriers, not tissue erasers.  First assign all
         * non-guide pixels from their topology component, then give pixels under
         * the guide stroke to the nearest side.  This keeps the final partition
         * exhaustive while preserving the guide-defined boundary. */
        int[] pixelOwner = new int[outputMask.length];
        Arrays.fill(pixelOwner, -1);
        for (int pixel = 0; pixel < outputMask.length; pixel++) {
            if ((outputMask[pixel] & 0xff) == 0) continue;
            int topologyId = labeling.labels[pixel];
            if (topologyId > 0 && topologyId < ownerByTopologyId.length)
                pixelOwner[pixel] = ownerByTopologyId[topologyId];
        }
        restoreGuidePixels(outputMask, pixelOwner, width, height,
            Math.max(4, barrierWidthPixels + 3));

        int[] piecesPerSeed = new int[seeds.size()];
        long acceptedPixels = 0;
        for (int seedIndex = 0; seedIndex < seeds.size(); seedIndex++) {
            byte[] oneOwner = new byte[outputMask.length];
            for (int pixel = 0; pixel < oneOwner.length; pixel++) {
                if (pixelOwner[pixel] == seedIndex) oneOwner[pixel] = (byte) 255;
            }
            Labeling finalPieces = labelComponents(oneOwner, null, width, height,
                dynamicMinimum, null, crop);
            int largestPiece = 0;
            for (Component component : finalPieces.components)
                if (component.area > largestPiece) largestPiece = component.area;
            int meaningfulPieceMinimum = Math.max(dynamicMinimum,
                (int) Math.round(largestPiece * 0.02));
            Collections.sort(finalPieces.components, new Comparator<Component>() {
                public int compare(Component left, Component right) {
                    int byY = Double.compare(left.centerY(), right.centerY());
                    return byY != 0 ? byY : Double.compare(left.centerX(), right.centerX());
                }
            });
            int part = 1;
            long filteredPixels = 0;
            int filteredPieces = 0;
            for (Component component : finalPieces.components) {
                if (component.area < meaningfulPieceMinimum) {
                    filteredPixels += component.area;
                    filteredPieces++;
                    continue;
                }
                Roi roi;
                if ("PL".equals(type)) {
                    roi = componentToCenterline(finalPieces.labels, component, oneOwner, width, crop);
                    if (roi == null) {
                        result.warnings.add("PL output component " + component.id +
                            " for " + seeds.get(seedIndex).region +
                            " could not be converted to a centerline and was omitted.");
                        continue;
                    }
                } else {
                    roi = componentToRoi(finalPieces.labels, component, oneOwner, width, crop);
                }
                if (roi == null) continue;
                roi.setName(AUTO_PREFIX + type + "_" + seeds.get(seedIndex).region +
                    "_part" + part++);
                applyPreviewStyle(roi, type);
                result.candidates.add(roi);
                piecesPerSeed[seedIndex]++;
                acceptedPixels += component.area;
            }
            if (filteredPieces > 0)
                result.warnings.add(type + " " + seeds.get(seedIndex).region +
                    ": filtered " + filteredPieces + " detached fragment(s), " +
                    filteredPixels + " pixels total (<2% of the largest part). " +
                    "Lower Minimum component area only if those fragments are true target tissue.");
        }
        for (int s = 0; s < seeds.size(); s++) {
            if (piecesPerSeed[s] == 0)
                result.warnings.add(type + " has no candidate for " + seeds.get(s).region +
                    ". Add/extend a fissure or move the seed.");
        }
        return acceptedPixels;
    }

    private static void restoreGuidePixels(byte[] outputMask, int[] pixelOwner,
            int width, int height, int maximumRadius) {
        int[] originalOwners = pixelOwner.clone();
        for (int pixel = 0; pixel < outputMask.length; pixel++) {
            if ((outputMask[pixel] & 0xff) == 0 || originalOwners[pixel] >= 0) continue;
            int x = pixel % width;
            int y = pixel / width;
            int owner = nearestOwnerOnRing(originalOwners, x, y, width, height, maximumRadius);
            if (owner >= 0) pixelOwner[pixel] = owner;
        }
    }

    private static int nearestOwnerOnRing(int[] owners, int x, int y,
            int width, int height, int maximumRadius) {
        for (int radius = 1; radius <= maximumRadius; radius++) {
            int bestOwner = Integer.MAX_VALUE;
            for (int yy = y - radius; yy <= y + radius; yy++) {
                for (int xx = x - radius; xx <= x + radius; xx++) {
                    if (xx < 0 || yy < 0 || xx >= width || yy >= height) continue;
                    if (Math.abs(xx - x) != radius && Math.abs(yy - y) != radius) continue;
                    int owner = owners[yy * width + xx];
                    if (owner >= 0 && owner < bestOwner) bestOwner = owner;
                }
            }
            if (bestOwner != Integer.MAX_VALUE) return bestOwner;
        }
        return -1;
    }

    private static int nearestAllowedSeed(Component component, Labeling layer, List<Seed> seeds,
            Labeling parentLabels, int[] seedParents, int width, int height) {
        int parent = componentParent(component, layer, parentLabels, width, height);
        int bestSeed = -1;
        double best = Double.POSITIVE_INFINITY;
        for (int s = 0; s < seeds.size(); s++) {
            if (parent >= 0 && seedParents[s] >= 0 && parent != seedParents[s]) continue;
            if (component.minSeedDistanceSquared[s] < best) {
                best = component.minSeedDistanceSquared[s];
                bestSeed = s;
            }
        }
        return bestSeed;
    }

    private static long countComponentPixels(int[] labels, Component component, byte[] outputMask,
            int fullWidth) {
        long count = 0;
        for (int y = component.minY; y <= component.maxY; y++) {
            int row = y * fullWidth;
            for (int x = component.minX; x <= component.maxX; x++) {
                int pixel = row + x;
                if (labels[pixel] == component.id && (outputMask[pixel] & 0xff) != 0) count++;
            }
        }
        return count;
    }

    private static Roi componentToRoi(int[] labels, Component component, byte[] outputMask,
            int fullWidth, Rectangle crop) {
        int width = component.maxX - component.minX + 1;
        int height = component.maxY - component.minY + 1;
        if (width <= 0 || height <= 0) return null;
        byte[] pixels = new byte[width * height];
        for (int y = component.minY; y <= component.maxY; y++) {
            int sourceRow = y * fullWidth;
            int targetRow = (y - component.minY) * width;
            for (int x = component.minX; x <= component.maxX; x++) {
                int source = sourceRow + x;
                if (labels[source] == component.id && (outputMask[source] & 0xff) != 0)
                    pixels[targetRow + x - component.minX] = (byte) 255;
            }
        }
        ByteProcessor processor = new ByteProcessor(width, height, pixels, null);
        processor.setThreshold(1, 255, ImageProcessor.NO_LUT_UPDATE);
        Roi roi = new ThresholdToSelection().convert(processor);
        if (roi == null) return null;
        Rectangle bounds = roi.getBounds();
        roi.setLocation(bounds.x + crop.x + component.minX, bounds.y + crop.y + component.minY);
        return roi;
    }

    /** Convert a narrow PL band into one polyline by skeletonizing it and taking
     * the longest geodesic path.  This preserves legacy CRA's PL Length semantics
     * instead of silently exporting a traced area ROI. */
    private static Roi componentToCenterline(int[] labels, Component component,
            byte[] outputMask, int fullWidth, Rectangle crop) {
        int width = component.maxX - component.minX + 3;
        int height = component.maxY - component.minY + 3;
        if (width <= 2 || height <= 2) return null;
        byte[] band = new byte[width * height];
        for (int y = component.minY; y <= component.maxY; y++) {
            int sourceRow = y * fullWidth;
            int localRow = (y - component.minY + 1) * width;
            for (int x = component.minX; x <= component.maxX; x++) {
                int source = sourceRow + x;
                if (labels[source] == component.id && (outputMask[source] & 0xff) != 0)
                    band[localRow + x - component.minX + 1] = (byte) 255;
            }
        }
        ByteProcessor source = new ByteProcessor(width, height, band, null);
        source.invert(); // BinaryProcessor skeletonizes black foreground on white background.
        BinaryProcessor skeleton = new BinaryProcessor(source);
        skeleton.skeletonize();
        byte[] pixels = (byte[]) skeleton.getPixels();

        int first = -1;
        int firstEndpoint = -1;
        int skeletonPixels = 0;
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                int index = y * width + x;
                if ((pixels[index] & 0xff) != 0) continue;
                skeletonPixels++;
                if (first < 0) first = index;
                if (skeletonDegree(pixels, index, width, height) == 1 && firstEndpoint < 0)
                    firstEndpoint = index;
            }
        }
        int start = firstEndpoint >= 0 ? firstEndpoint : first;
        if (start < 0) return null;
        DijkstraPath firstPass = skeletonDistances(pixels, width, height, start, false);
        DijkstraPath secondPass = skeletonDistances(pixels, width, height, firstPass.farthest, true);
        if (secondPass.farthest < 0) return null;

        List<Integer> reversePath = new ArrayList<Integer>();
        int cursor = secondPass.farthest;
        while (cursor >= 0) {
            reversePath.add(Integer.valueOf(cursor));
            if (cursor == firstPass.farthest) break;
            cursor = secondPass.predecessor[cursor];
        }
        if (reversePath.size() < 2) return null;
        /* A strongly branched band would lose real PL length if only its longest
         * path were kept.  Omit it rather than return a plausible but incomplete
         * line; the missing-candidate gate then forces manual review. */
        if (skeletonPixels > 0 && reversePath.size() / (double) skeletonPixels < 0.85)
            return null;
        Collections.reverse(reversePath);
        /* Pixel skeletons contain staircase zig-zags; moderate geometric
         * simplification recovers the smooth traced-line length more faithfully
         * than exporting every raster step. */
        List<Integer> simplified = simplifyPixelPath(reversePath, width, 1.5);
        if (simplified.size() < 2) return null;
        float[] xs = new float[simplified.size()];
        float[] ys = new float[simplified.size()];
        for (int i = 0; i < simplified.size(); i++) {
            int pixel = simplified.get(i).intValue();
            int x = pixel % width;
            int y = pixel / width;
            xs[i] = x - 1 + component.minX + crop.x;
            ys[i] = y - 1 + component.minY + crop.y;
        }
        return new PolygonRoi(xs, ys, xs.length, Roi.POLYLINE);
    }

    private static int skeletonDegree(byte[] pixels, int index, int width, int height) {
        int x = index % width;
        int y = index / width;
        int degree = 0;
        for (int yy = y - 1; yy <= y + 1; yy++) {
            for (int xx = x - 1; xx <= x + 1; xx++) {
                if (xx == x && yy == y) continue;
                if (xx < 0 || yy < 0 || xx >= width || yy >= height) continue;
                if ((pixels[yy * width + xx] & 0xff) == 0) degree++;
            }
        }
        return degree;
    }

    private static final class DijkstraPath {
        int farthest = -1;
        int[] predecessor;
    }

    private static DijkstraPath skeletonDistances(byte[] pixels, int width, int height,
            int start, boolean keepPredecessor) {
        float[] distance = new float[pixels.length];
        Arrays.fill(distance, Float.POSITIVE_INFINITY);
        int[] predecessor = new int[pixels.length];
        Arrays.fill(predecessor, -1);
        MinHeap heap = new MinHeap(1024);
        distance[start] = 0f;
        heap.push(start, 0f);
        final int[] dx = {-1, 0, 1, -1, 1, -1, 0, 1};
        final int[] dy = {-1, -1, -1, 0, 0, 1, 1, 1};
        final float diagonal = 1.41421356f;
        int farthest = start;
        while (!heap.isEmpty()) {
            int current = heap.popIndex();
            float currentDistance = heap.lastPoppedKey();
            if (currentDistance != distance[current]) continue;
            if (currentDistance > distance[farthest]) farthest = current;
            int x = current % width;
            int y = current / width;
            for (int direction = 0; direction < 8; direction++) {
                int xx = x + dx[direction];
                int yy = y + dy[direction];
                if (xx < 0 || yy < 0 || xx >= width || yy >= height) continue;
                int next = yy * width + xx;
                if ((pixels[next] & 0xff) != 0) continue;
                float step = dx[direction] == 0 || dy[direction] == 0 ? 1f : diagonal;
                float candidate = currentDistance + step;
                if (candidate < distance[next]) {
                    distance[next] = candidate;
                    predecessor[next] = current;
                    heap.push(next, candidate);
                }
            }
        }
        DijkstraPath answer = new DijkstraPath();
        answer.farthest = farthest;
        answer.predecessor = keepPredecessor ? predecessor : null;
        return answer;
    }

    private static List<Integer> simplifyPixelPath(List<Integer> path, int width, double tolerance) {
        if (path.size() <= 2) return path;
        boolean[] keep = new boolean[path.size()];
        keep[0] = true;
        keep[path.size() - 1] = true;
        int[] rangeStart = new int[path.size()];
        int[] rangeEnd = new int[path.size()];
        int stack = 0;
        rangeStart[stack] = 0;
        rangeEnd[stack++] = path.size() - 1;
        double toleranceSquared = tolerance * tolerance;
        while (stack > 0) {
            int first = rangeStart[--stack];
            int last = rangeEnd[stack];
            int a = path.get(first).intValue();
            int b = path.get(last).intValue();
            double ax = a % width;
            double ay = a / width;
            double bx = b % width;
            double by = b / width;
            double best = -1.0;
            int bestIndex = -1;
            for (int i = first + 1; i < last; i++) {
                int p = path.get(i).intValue();
                double distance = pointSegmentDistanceSquared(p % width, p / width, ax, ay, bx, by);
                if (distance > best) {
                    best = distance;
                    bestIndex = i;
                }
            }
            if (bestIndex >= 0 && best > toleranceSquared) {
                keep[bestIndex] = true;
                rangeStart[stack] = first;
                rangeEnd[stack++] = bestIndex;
                rangeStart[stack] = bestIndex;
                rangeEnd[stack++] = last;
            }
        }
        List<Integer> simplified = new ArrayList<Integer>();
        for (int i = 0; i < path.size(); i++) if (keep[i]) simplified.add(path.get(i));
        return simplified;
    }

    private static double pointSegmentDistanceSquared(double px, double py,
            double ax, double ay, double bx, double by) {
        double dx = bx - ax;
        double dy = by - ay;
        if (dx == 0.0 && dy == 0.0) {
            dx = px - ax;
            dy = py - ay;
            return dx * dx + dy * dy;
        }
        double t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy);
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        double qx = ax + t * dx;
        double qy = ay + t * dy;
        dx = px - qx;
        dy = py - qy;
        return dx * dx + dy * dy;
    }

    private interface PixelRule {
        boolean include(int pixelIndex);
    }

    private static void addDerivedTotal(String canonicalName, int width, int height,
            Rectangle crop, PixelRule rule, List<Roi> output) {
        byte[] mask = new byte[width * height];
        for (int p = 0; p < mask.length; p++) {
            if (rule.include(p)) mask[p] = (byte) 255;
        }
        int minimum = (int) Math.max(1, Math.round(countMask(mask) * 0.0001));
        mask = removeSmallComponents(mask, width, height, crop, minimum);
        Roi roi = maskToRoi(mask, width, height, crop);
        if (roi == null) return;
        roi.setName(AUTO_PREFIX + canonicalName);
        roi.setStrokeColor(PREVIEW_TOTAL_COLOR);
        roi.setStrokeWidth(1.5);
        output.add(roi);
    }

    private static byte[] removeSmallComponents(byte[] mask, int width, int height,
            Rectangle crop, int minimumArea) {
        Labeling labeling = labelComponents(mask, null, width, height, minimumArea, null, crop);
        int largestId = 0;
        for (Component component : labeling.components) if (component.id > largestId) largestId = component.id;
        boolean[] keep = new boolean[largestId + 1];
        for (Component component : labeling.components) keep[component.id] = true;
        byte[] filtered = new byte[mask.length];
        for (int i = 0; i < filtered.length; i++) {
            int id = labeling.labels[i];
            if (id > 0 && id < keep.length && keep[id]) filtered[i] = (byte) 255;
        }
        return filtered;
    }

    private static void createAreaCandidates(String type, List<Seed> seeds, byte[] owner,
            Rectangle crop, int minimumArea, PixelRule rule, List<Roi> output) {
        int width = crop.width;
        int height = crop.height;
        byte[][] masks = new byte[seeds.size()][width * height];
        for (int p = 0; p < owner.length; p++) {
            int encodedOwner = owner[p] & 0xff;
            if (encodedOwner == 0 || !rule.include(p)) continue;
            masks[encodedOwner - 1][p] = (byte) 255;
        }

        for (int seedIndex = 0; seedIndex < seeds.size(); seedIndex++) {
            Roi selection = maskToRoi(masks[seedIndex], width, height, crop);
            if (selection == null) continue;
            List<Roi> pieces = splitComponents(selection, minimumArea);
            Collections.sort(pieces, new Comparator<Roi>() {
                public int compare(Roi left, Roi right) {
                    double[] lc = left.getContourCentroid();
                    double[] rc = right.getContourCentroid();
                    int byY = Double.compare(lc[1], rc[1]);
                    return byY != 0 ? byY : Double.compare(lc[0], rc[0]);
                }
            });
            for (int part = 0; part < pieces.size(); part++) {
                Roi candidate = pieces.get(part);
                candidate.setName(AUTO_PREFIX + type + "_" + seeds.get(seedIndex).region + "_part" + (part + 1));
                applyPreviewStyle(candidate, type);
                output.add(candidate);
            }
        }
    }

    private static Roi maskToRoi(byte[] mask, int width, int height, Rectangle crop) {
        int minX = width, minY = height, maxX = -1, maxY = -1;
        for (int y = 0; y < height; y++) {
            int row = y * width;
            for (int x = 0; x < width; x++) {
                if ((mask[row + x] & 0xff) == 0) continue;
                if (x < minX) minX = x;
                if (x > maxX) maxX = x;
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
            }
        }
        if (maxX < minX || maxY < minY) return null;
        int subWidth = maxX - minX + 1;
        int subHeight = maxY - minY + 1;
        byte[] sub = new byte[subWidth * subHeight];
        for (int y = 0; y < subHeight; y++)
            System.arraycopy(mask, (minY + y) * width + minX, sub, y * subWidth, subWidth);
        ByteProcessor processor = new ByteProcessor(subWidth, subHeight, sub, null);
        processor.setThreshold(1, 255, ImageProcessor.NO_LUT_UPDATE);
        Roi roi = new ThresholdToSelection().convert(processor);
        if (roi == null) return null;
        Rectangle bounds = roi.getBounds();
        roi.setLocation(bounds.x + crop.x + minX, bounds.y + crop.y + minY);
        return roi;
    }

    private static List<Roi> splitComponents(Roi selection, int minimumArea) {
        Roi[] raw = selection instanceof ShapeRoi ? ((ShapeRoi) selection).getRois() : new Roi[]{selection};
        List<Roi> pieces = new ArrayList<Roi>();
        for (Roi piece : raw) {
            if (!piece.isArea()) continue;
            if (piece.getStatistics().area < minimumArea) continue;
            pieces.add((Roi) piece.clone());
        }
        if (pieces.isEmpty() && selection.isArea() && selection.getStatistics().area >= minimumArea)
            pieces.add((Roi) selection.clone());
        return pieces;
    }

    /** Split PL into per-lobule candidates by walking PL_TOTAL's ORIGINAL traced
     * vertices and cutting wherever pixel ownership changes, instead of rasterizing
     * and skeletonizing.  This preserves the exact hand-traced coordinates (true
     * PL length, no raster/skeleton round trip) and only requires PL_TOTAL to be a
     * line; it does not need per-lobule closed loops.  Ownership reuses the same
     * barrier + nearest-seed component labeling as the ML/GL split, so PL, ML and
     * GL boundaries stay consistent with each other. */
    private static long splitPLDirect(Roi plTotal, byte[] plMask, byte[] barriers,
            List<Seed> seeds, Labeling parentLabels, int[] seedParents, Rectangle crop,
            int requestedMinimumArea, int barrierWidthPixels, Result result) {
        int width = crop.width;
        int height = crop.height;
        long layerArea = countMask(plMask);
        int dynamicMinimum = (int) Math.max(requestedMinimumArea, Math.round(layerArea * 0.0001));
        Labeling labeling = labelComponents(plMask, barriers, width, height, dynamicMinimum, seeds, crop);
        int maximumTopologyId = 0;
        for (Component component : labeling.components)
            if (component.id > maximumTopologyId) maximumTopologyId = component.id;
        int[] ownerByTopologyId = new int[maximumTopologyId + 1];
        Arrays.fill(ownerByTopologyId, -1);

        for (int c = 0; c < labeling.components.size(); c++) {
            Component component = labeling.components.get(c);
            int parent = componentParent(component, labeling, parentLabels, width, height);
            int bestSeed = -1;
            double best = Double.POSITIVE_INFINITY;
            double second = Double.POSITIVE_INFINITY;
            for (int s = 0; s < seeds.size(); s++) {
                if (parent >= 0 && seedParents[s] >= 0 && parent != seedParents[s]) continue;
                double distance = component.minSeedDistanceSquared[s];
                if (distance < best) {
                    second = best;
                    best = distance;
                    bestSeed = s;
                } else if (distance < second) {
                    second = distance;
                }
            }
            if (bestSeed < 0) {
                result.warnings.add("PL component " + component.id +
                    " could not be assigned within its tissue parent.");
                continue;
            }
            if (component.id < ownerByTopologyId.length)
                ownerByTopologyId[component.id] = bestSeed;
            result.diagnostics.add("PL component " + component.id + ": area=" + component.area +
                ", parent=" + parent + ", owner=" + seeds.get(bestSeed).region +
                ", nearestDistance=" + String.format("%.1f px", Math.sqrt(best)));
            if (second < Double.POSITIVE_INFINITY) {
                double margin = Math.sqrt(second) - Math.sqrt(best);
                if (margin < 20.0)
                    result.warnings.add("PL component " + component.id +
                        " has a low seed-distance margin (" + String.format("%.1f px", margin) +
                        "); inspect its boundary carefully.");
            }
        }

        int[] pixelOwner = new int[plMask.length];
        Arrays.fill(pixelOwner, -1);
        for (int pixel = 0; pixel < plMask.length; pixel++) {
            if ((plMask[pixel] & 0xff) == 0) continue;
            int topologyId = labeling.labels[pixel];
            if (topologyId > 0 && topologyId < ownerByTopologyId.length)
                pixelOwner[pixel] = ownerByTopologyId[topologyId];
        }
        restoreGuidePixels(plMask, pixelOwner, width, height, Math.max(4, barrierWidthPixels + 3));

        long assignedPixels = 0;
        for (int pixel = 0; pixel < plMask.length; pixel++) {
            if ((plMask[pixel] & 0xff) != 0 && pixelOwner[pixel] >= 0) assignedPixels++;
        }

        emitPLVertexSegments(plTotal, pixelOwner, width, height, crop, seeds, result);
        return assignedPixels;
    }

    private static void emitPLVertexSegments(Roi plTotal, int[] pixelOwner, int width, int height,
            Rectangle crop, List<Seed> seeds, Result result) {
        FloatPolygon points = plTotal.getFloatPolygon();
        if (points == null || points.npoints < 2) {
            result.warnings.add("PL_TOTAL has no usable vertices; PL candidates were not generated.");
            return;
        }
        int[] partNumbers = new int[seeds.size()];
        Arrays.fill(partNumbers, 1);
        int[] piecesPerSeed = new int[seeds.size()];

        int start = 0;
        int previousOwner = ownerAt(points.xpoints[0], points.ypoints[0], pixelOwner, width, height, crop);
        for (int i = 1; i <= points.npoints; i++) {
            int owner = i < points.npoints
                ? ownerAt(points.xpoints[i], points.ypoints[i], pixelOwner, width, height, crop)
                : -2;
            if (owner != previousOwner) {
                if (previousOwner >= 0 && i - start >= 2) {
                    int count = i - start;
                    float[] xs = new float[count];
                    float[] ys = new float[count];
                    System.arraycopy(points.xpoints, start, xs, 0, count);
                    System.arraycopy(points.ypoints, start, ys, 0, count);
                    double segmentLength = polylineLength(xs, ys);
                    if (segmentLength < 3.0) {
                        result.diagnostics.add("PL micro-segment for " + seeds.get(previousOwner).region +
                            " skipped (" + String.format("%.2f px", segmentLength) + ", below 3 px noise threshold).");
                    } else {
                        PolygonRoi line = new PolygonRoi(xs, ys, count, Roi.POLYLINE);
                        line.setName(AUTO_PREFIX + "PL_" + seeds.get(previousOwner).region +
                            "_part" + partNumbers[previousOwner]++);
                        applyPreviewStyle(line, "PL");
                        result.candidates.add(line);
                        piecesPerSeed[previousOwner]++;
                    }
                }
                start = i;
                previousOwner = owner;
            }
        }
        for (int s = 0; s < seeds.size(); s++) {
            if (piecesPerSeed[s] == 0)
                result.warnings.add("PL has no candidate for " + seeds.get(s).region +
                    ". Add/extend a fissure or move the seed.");
        }
    }

    private static int ownerAt(float x, float y, int[] pixelOwner, int width, int height, Rectangle crop) {
        int gx = (int) Math.round(x - crop.x);
        int gy = (int) Math.round(y - crop.y);
        gx = Math.max(0, Math.min(width - 1, gx));
        gy = Math.max(0, Math.min(height - 1, gy));
        int direct = pixelOwner[gy * width + gx];
        if (direct >= 0) return direct;
        return nearestOwnerOnRing(pixelOwner, gx, gy, width, height, 8);
    }

    private static double polylineLength(float[] xs, float[] ys) {
        double length = 0.0;
        for (int i = 1; i < xs.length; i++) {
            double dx = xs[i] - xs[i - 1];
            double dy = ys[i] - ys[i - 1];
            length += Math.sqrt(dx * dx + dy * dy);
        }
        return length;
    }

    public static Map<String, Roi> toNamedMap(Roi[] rois) {
        Map<String, Roi> result = new LinkedHashMap<String, Roi>();
        for (Roi roi : rois) {
            if (roi != null && roi.getName() != null) result.put(roi.getName(), roi);
        }
        return result;
    }

    private static final class MinHeap {
        private int[] indexes;
        private float[] keys;
        private int size;
        private float poppedKey;

        MinHeap(int initialCapacity) {
            indexes = new int[initialCapacity];
            keys = new float[initialCapacity];
        }

        boolean isEmpty() { return size == 0; }
        float lastPoppedKey() { return poppedKey; }

        void push(int index, float key) {
            ensureCapacity(size + 1);
            int position = size++;
            while (position > 0) {
                int parent = (position - 1) >>> 1;
                if (keys[parent] <= key) break;
                indexes[position] = indexes[parent];
                keys[position] = keys[parent];
                position = parent;
            }
            indexes[position] = index;
            keys[position] = key;
        }

        int popIndex() {
            int answer = indexes[0];
            poppedKey = keys[0];
            int lastIndex = indexes[--size];
            float lastKey = keys[size];
            if (size > 0) {
                int position = 0;
                while (true) {
                    int left = position * 2 + 1;
                    if (left >= size) break;
                    int right = left + 1;
                    int child = right < size && keys[right] < keys[left] ? right : left;
                    if (keys[child] >= lastKey) break;
                    indexes[position] = indexes[child];
                    keys[position] = keys[child];
                    position = child;
                }
                indexes[position] = lastIndex;
                keys[position] = lastKey;
            }
            return answer;
        }

        private void ensureCapacity(int needed) {
            if (needed <= indexes.length) return;
            int next = Math.max(needed, indexes.length * 2);
            indexes = Arrays.copyOf(indexes, next);
            keys = Arrays.copyOf(keys, next);
        }
    }
}
