package cra.autosplit;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import ij.process.FloatPolygon;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

/** ImageJ 1 user interface and ROI Manager integration for CRA AutoSplit. */
public final class CRAAutoSplitPlugin implements PlugIn {
    private static final String TITLE = "CRA AutoSplit A003 alpha4";
    private static final String SEED_PREFIX = "SEED_";
    private static final String DIVIDER_PREFIX = "DIVIDER_";

    private static final String[] REGION_LABELS = {
        "1Cb", "2Cb", "3Cb", "4/5Cb", "6Cb", "7Cb", "8Cb", "9Cb", "10Cb",
        "Sim", "Crus1", "Crus2", "PM", "Cop", "PFl", "Fl", "Parafl", "Other..."
    };

    private static String lastRegionLabel = "Crus1";
    private static int nextDividerNumber = 1;
    private static int optionBarrierWidth = 7;
    private static int optionMinimumComponent = 500;
    /* When PL_TOTAL is a line, PL candidates are cut directly from its ORIGINAL
     * traced vertices (no rasterize/skeletonize round trip), using the same
     * barrier + nearest-seed ownership as ML/GL, so PL length is preserved
     * exactly.  When PL_TOTAL is an area/band, the lossier skeleton-centreline
     * fallback is used instead and stays labeled EXPERIMENTAL in the dialog.
     * Off by default: opt in once the SEED_-anchored workflow (not local-PL
     * derived anchors, which this option is blocked for) has been verified. */
    private static boolean optionCreatePL = false;
    private static boolean optionCreateML = true;
    private static boolean optionCreateGL = true;
    private static boolean optionCreateTotals = true;
    private static final Set<String> lastPreviewNames = new HashSet<String>();
    private static final Set<String> lastPreviewRegions = new HashSet<String>();
    private static double lastPreviewCoverage = -1.0;
    private static boolean lastPreviewPL;
    private static boolean lastPreviewML;
    private static boolean lastPreviewGL;

    public void run(String argument) {
        String command = argument == null ? "guide" : argument.trim().toLowerCase(Locale.US);
        try {
            if ("start".equals(command) || "smartstart".equals(command)) {
                smartStart();
            } else if ("addseedset".equals(command)) {
                addSeedSet();
            } else if ("addseed".equals(command)) {
                addSeed();
            } else if ("adddivider".equals(command)) {
                addDivider();
            } else if ("validate".equals(command) || "validateseeds".equals(command)) {
                validate();
            } else if ("preview".equals(command)) {
                preview();
            } else if ("accept".equals(command)) {
                acceptPreview();
            } else if ("clear".equals(command) || "clearpreview".equals(command)) {
                clearPreview(true);
            } else if ("deletehelpers".equals(command)) {
                deleteHelpers(true);
            } else {
                showGuide();
            }
        } catch (IllegalArgumentException expected) {
            IJ.error(TITLE, expected.getMessage());
        } catch (Throwable unexpected) {
            IJ.log(TITLE + " unexpected error: " + unexpected);
            IJ.handleException(unexpected);
        }
    }

    /** One entry point for routine use: inspect first, then open Preview only
     * when every blocking requirement is satisfied. */
    private static void smartStart() {
        ImagePlus image = requireImage();
        RoiManager manager = RoiManager.getRoiManager();
        Validation validation = inspect(manager.getRoisAsArray());
        if (!validation.errors.isEmpty()) {
            IJ.showMessage(TITLE + " - Input check", validationReport(validation, false));
            image.updateAndDraw();
            return;
        }
        preview();
    }

    private static void addSeedSet() {
        ImagePlus image = requireImage();
        Roi active = image.getRoi();
        if (!(active instanceof PointRoi))
            throw new IllegalArgumentException(
                "Multi-point 도구로 보이는 lobule을 해부학적 순서대로 한 번씩 클릭한 뒤 실행하세요.");
        FloatPolygon points = active.getFloatPolygon();
        if (points == null || points.npoints < 2)
            throw new IllegalArgumentException("Seed set에는 최소 두 개의 Point가 필요합니다.");

        GenericDialog dialog = new GenericDialog(TITLE + " - Add Seed Set");
        dialog.addMessage("클릭한 순서와 정확히 같은 순서로 lobule 이름을 쉼표로 입력합니다.\n" +
            "예: Sim,Crus1,Crus2,PM,PFl,Fl\n" +
            "저장 시 기존 SEED_ helpers는 이 새 set으로 교체됩니다.");
        dialog.addStringField("Lobules in click order", "Sim,Crus1,Crus2,PM,PFl,Fl", 42);
        dialog.setOKLabel("Create seed set");
        dialog.showDialog();
        if (dialog.wasCanceled()) return;
        String[] raw = dialog.getNextString().split("[,;]");
        if (raw.length != points.npoints)
            throw new IllegalArgumentException("Point count (" + points.npoints +
                ") and lobule-name count (" + raw.length + ") do not match.");
        List<String> regions = new ArrayList<String>();
        Set<String> unique = new HashSet<String>();
        for (String value : raw) {
            String region = sanitizeRegion(value);
            if (!unique.add(region))
                throw new IllegalArgumentException("Duplicate lobule name in seed set: " + region);
            regions.add(region);
        }

        RoiManager manager = RoiManager.getRoiManager();
        deleteByPrefixes(manager, new String[]{SEED_PREFIX});
        for (int i = 0; i < regions.size(); i++) {
            PointRoi seed = new PointRoi(points.xpoints[i], points.ypoints[i]);
            seed.setName(SEED_PREFIX + regions.get(i));
            seed.setPointType(PointRoi.CROSS);
            seed.setSize(3);
            seed.setStrokeColor(Color.CYAN);
            manager.addRoi(seed);
        }
        image.deleteRoi();
        showAllWithLabels(manager, image);
        IJ.showMessage(TITLE, regions.size() +
            " one-point seed helpers were created. Run Validate Inputs next.");
    }

    private static void addSeed() {
        ImagePlus image = requireImage();
        Roi active = image.getRoi();
        if (!(active instanceof PointRoi)) {
            throw new IllegalArgumentException(
                "Point tool로 해당 lobule 내부를 클릭한 뒤 다시 실행하세요.\n" +
                "여러 점이 보이더라도 이 명령은 마지막 점 하나만 안전하게 저장합니다.");
        }
        FloatPolygon points = active.getFloatPolygon();
        if (points == null || points.npoints < 1)
            throw new IllegalArgumentException("활성 Point ROI에서 좌표를 읽을 수 없습니다.");

        GenericDialog dialog = new GenericDialog(TITLE + " - Add Lobule Seed");
        dialog.addMessage("각 lobule마다 내부의 대표 위치 한 곳을 지정합니다.\n" +
            "기존 같은 이름의 seed는 이 점 하나로 교체됩니다.");
        dialog.addChoice("Lobule", REGION_LABELS, lastRegionLabel);
        dialog.addStringField("Other name (only for Other...)", "", 18);
        dialog.setOKLabel("Save one seed");
        dialog.showDialog();
        if (dialog.wasCanceled()) return;
        String choice = dialog.getNextChoice();
        String custom = dialog.getNextString();
        String region = "Other...".equals(choice) ? custom : choice;
        region = sanitizeRegion(region);
        lastRegionLabel = choice;

        int last = points.npoints - 1;
        PointRoi seed = new PointRoi(points.xpoints[last], points.ypoints[last]);
        seed.setName(SEED_PREFIX + region);
        seed.setPointType(PointRoi.CROSS);
        seed.setSize(3);
        seed.setStrokeColor(Color.CYAN);

        RoiManager manager = RoiManager.getRoiManager();
        deleteExactName(manager, seed.getName());
        manager.addRoi(seed);
        image.setRoi((Roi) seed.clone());
        showAllWithLabels(manager, image);
        IJ.showStatus(TITLE + ": saved " + seed.getName() + " at one point");
    }

    private static void addDivider() {
        ImagePlus image = requireImage();
        Roi active = image.getRoi();
        if (active == null || !active.isLine() || active instanceof PointRoi)
            throw new IllegalArgumentException(
                "Segmented Line 도구로 필요한 경계를 그린 뒤 다시 실행하세요.\n" +
                "Divider는 분리할 두 cortical region 사이를 완전히 가로질러야 합니다.");

        GenericDialog dialog = new GenericDialog(TITLE + " - Add Divider Guide");
        dialog.addMessage("Fissure만으로 분리가 안 되는 곳에만 사용합니다.\n" +
            "해부학적 측정 ROI가 아니라 자동 분할용 보조선입니다.");
        dialog.addStringField("Divider label", String.valueOf(nextDividerNumber), 14);
        dialog.setOKLabel("Add divider");
        dialog.showDialog();
        if (dialog.wasCanceled()) return;
        String token = sanitizeRegion(dialog.getNextString());
        String name = DIVIDER_PREFIX + token;
        nextDividerNumber++;

        Roi divider = (Roi) active.clone();
        divider.setName(name);
        divider.setStrokeColor(new Color(255, 153, 0));
        divider.setStrokeWidth(2.0);
        RoiManager manager = RoiManager.getRoiManager();
        deleteExactName(manager, name);
        manager.addRoi(divider);
        image.setRoi((Roi) divider.clone());
        showAllWithLabels(manager, image);
        IJ.showStatus(TITLE + ": added " + name);
    }

    private static void validate() {
        ImagePlus image = requireImage();
        RoiManager manager = RoiManager.getRoiManager();
        Validation validation = inspect(manager.getRoisAsArray());
        IJ.showMessage(TITLE + " - Validation", validationReport(validation, true));
        image.updateAndDraw();
    }

    private static String validationReport(Validation validation, boolean includeReadyHint) {
        StringBuilder report = new StringBuilder();
        report.append(validation.errors.isEmpty() ? "VALIDATION: READY\n" : "VALIDATION: NOT READY\n");
        report.append("Seeds: ").append(validation.seeds.size()).append('\n');
        report.append("Lobule anchor source: ").append(validation.anchorSource).append('\n');
        report.append("Fissure guides: ").append(validation.fissureCount).append('\n');
        report.append("Divider guides: ").append(validation.dividerCount).append('\n');
        appendMessages(report, "ERROR", validation.errors);
        appendMessages(report, "WARNING", validation.warnings);
        if (validation.errors.isEmpty() && includeReadyHint)
            report.append("\nPreview를 실행할 수 있습니다. 최종 해부학적 판단은 연구자가 합니다.");
        else if (validation.errors.isEmpty())
            report.append("\n입력 검사를 통과했습니다.");
        else
            report.append("\n위 ERROR만 수정하면 됩니다. 원본 ROI는 변경되지 않았습니다.");
        return report.toString();
    }

    private static void preview() {
        ImagePlus image = requireImage();
        RoiManager manager = RoiManager.getRoiManager();
        Validation validation = inspect(manager.getRoisAsArray());
        if (!validation.errors.isEmpty())
            throw new IllegalArgumentException(validationFailureMessage(validation));

        Roi plTotal = validation.namedRois.get("PL_TOTAL");
        boolean plTotalIsBand = plTotal != null && plTotal.isArea();
        if (validation.derivedFromLocalPL) optionCreatePL = false;

        GenericDialog dialog = new GenericDialog(TITLE + " - Generate Review Candidates");
        dialog.addMessage("원본 ROI는 수정하지 않습니다. 생성물에는 AUTO_REVIEW_ 접두사가 붙습니다.\n" +
            "표시색: PL=적색, ML=녹색, GL=황색, 파생 TOTAL=하늘색.\n" +
            "AUTO_REVIEW_는 임시 안전표시이며 Accept 때 일괄 제거됩니다.\n" +
            "Lobule anchor: " + validation.anchorSource + "\n" +
            (validation.derivedFromLocalPL
                ? "기존 local PL을 보존하고 ML/GL만 자동 생성하는 방식이 권장됩니다.\n"
                : "SEED_는 승인 시 자동 삭제할 수 있는 임시 이름 기준점입니다.\n") +
            (plTotalIsBand
                ? "주의: PL_TOTAL은 면적 band입니다. PL 중심선 생성은 실험 기능이며 기본 해제됩니다.\n"
                : "") +
            "모든 경계를 눈으로 확인한 뒤에만 Accept를 실행하세요.");
        dialog.addNumericField("Barrier width (pixels)", optionBarrierWidth, 0, 6, "fissure/divider blocking width");
        dialog.addNumericField("Minimum component area (pixels)", optionMinimumComponent, 0, 8, "remove tiny fragments");
        dialog.addCheckbox(plTotalIsBand
            ? "EXPERIMENTAL: derive PL centreline candidates from PL_TOTAL band"
            : "Create PL candidates (direct split of PL_TOTAL's traced vertices)", optionCreatePL);
        dialog.addCheckbox("Create ML candidates", optionCreateML);
        dialog.addCheckbox("Create GL candidates", optionCreateGL);
        dialog.addCheckbox("Derive missing ML_TOTAL / GL_TOTAL candidates", optionCreateTotals);
        dialog.setOKLabel("Generate preview");
        dialog.showDialog();
        if (dialog.wasCanceled()) return;

        optionBarrierWidth = (int) Math.round(dialog.getNextNumber());
        optionMinimumComponent = (int) Math.round(dialog.getNextNumber());
        optionCreatePL = dialog.getNextBoolean();
        optionCreateML = dialog.getNextBoolean();
        optionCreateGL = dialog.getNextBoolean();
        optionCreateTotals = dialog.getNextBoolean();
        if (dialog.invalidNumber())
            throw new IllegalArgumentException("Preview settings contain an invalid number.");
        if (validation.derivedFromLocalPL && optionCreatePL)
            throw new IllegalArgumentException(
                "Existing final local PL ROIs are being used as anchors.\n" +
                "Keep experimental PL generation OFF so those verified PL lines are preserved.\n" +
                "The plugin will generate ML/GL candidates only.");

        AutoSplitEngine.Options options = new AutoSplitEngine.Options();
        options.barrierWidthPixels = optionBarrierWidth;
        options.minimumComponentPixels = optionMinimumComponent;
        options.createPL = optionCreatePL;
        options.createML = optionCreateML;
        options.createGL = optionCreateGL;
        options.createDerivedTotals = optionCreateTotals;

        IJ.showStatus(TITLE + ": splitting PL, ML and GL independently...");
        IJ.showProgress(0.05);
        AutoSplitEngine.Result result = AutoSplitEngine.generate(
            cloneNamedRois(validation.namedRois), validation.seeds, options);
        IJ.showProgress(0.9);

        clearPreview(false);
        for (Roi candidate : result.candidates)
            manager.addRoi((Roi) candidate.clone());
        lastPreviewNames.clear();
        for (Roi candidate : result.candidates)
            lastPreviewNames.add(safeName(candidate));
        lastPreviewRegions.clear();
        for (AutoSplitEngine.Seed seed : validation.seeds)
            lastPreviewRegions.add(seed.region);
        lastPreviewCoverage = result.assignedPercent();
        lastPreviewPL = optionCreatePL;
        lastPreviewML = optionCreateML;
        lastPreviewGL = optionCreateGL;
        showAllWithLabels(manager, image);
        IJ.showProgress(1.0);

        StringBuilder message = new StringBuilder();
        message.append("Review candidates: ").append(result.candidates.size()).append('\n');
        message.append("Seeds: ").append(validation.seeds.size()).append('\n');
        message.append("Lobule anchor source: ").append(validation.anchorSource).append('\n');
        message.append("Cortical domain assigned: ")
            .append(String.format(Locale.US, "%.2f%%", result.assignedPercent())).append('\n');
        message.append("Layer overlap audit (pixels): PL-ML=").append(result.plMlOverlapPixels)
            .append(", PL-GL=").append(result.plGlOverlapPixels)
            .append(", ML-GL=").append(result.mlGlOverlapPixels).append('\n');
        message.append("Input nesting audit (pixels): GL+WM outside Cb=")
            .append(result.glwmOutsideCbPixels)
            .append(", WM outside GL+WM=").append(result.wmOutsideGlwmPixels)
            .append(", PL outside Cb=").append(result.plOutsideCbPixels)
            .append(", PL-WM overlap=").append(result.plWmOverlapPixels).append('\n');
        message.append("Runtime: ").append(result.elapsedMillis).append(" ms\n");
        appendMessages(message, "WARNING", result.warnings);
        message.append("\nAUTO_REVIEW_ boundaries are proposals, not final annotations.\n")
            .append("PL(red), ML(green), and GL(yellow) masks are mutually exclusive; adjacent outlines may touch.\n")
            .append("Accept removes AUTO_REVIEW_ from every approved ROI in one batch; manual renaming is unnecessary.\n")
            .append("Inspect every candidate. Add/extend DIVIDER_ if a boundary is wrong, then regenerate.");
        IJ.showMessage(TITLE + " - Preview ready", message.toString());
        IJ.showStatus(TITLE + ": preview ready; human review required");
    }

    private static void acceptPreview() {
        ImagePlus image = requireImage();
        RoiManager manager = RoiManager.getRoiManager();
        Roi[] rois = manager.getRoisAsArray();
        List<Integer> previewIndexes = new ArrayList<Integer>();
        Set<String> existing = new HashSet<String>();
        Set<String> future = new HashSet<String>();
        List<String> collisions = new ArrayList<String>();
        Validation validation = inspect(rois);
        List<String> seedRegions = new ArrayList<String>();
        for (AutoSplitEngine.Seed seed : validation.seeds)
            seedRegions.add(seed.region);

        for (int i = 0; i < rois.length; i++) {
            String name = safeName(rois[i]);
            if (name.startsWith(AutoSplitEngine.AUTO_PREFIX)) {
                previewIndexes.add(Integer.valueOf(i));
            } else if (!name.isEmpty()) {
                existing.add(name);
            }
        }
        if (previewIndexes.isEmpty())
            throw new IllegalArgumentException("No AUTO_REVIEW_ candidate exists. Run Preview first.");
        if (lastPreviewNames.isEmpty())
            throw new IllegalArgumentException(
                "The Preview manifest is not available (for example, Fiji was restarted).\n" +
                "Regenerate Preview before acceptance; this prevents stale or incomplete candidates from being finalized.");
        Set<String> currentPreviewNames = new HashSet<String>();
        for (Integer boxed : previewIndexes)
            currentPreviewNames.add(safeName(rois[boxed.intValue()]));
        if (!currentPreviewNames.equals(lastPreviewNames))
            throw new IllegalArgumentException(
                "AUTO_REVIEW candidates changed after Preview.\n" +
                "Clear and regenerate Preview before acceptance.");
        Set<String> currentRegions = new HashSet<String>();
        for (AutoSplitEngine.Seed seed : validation.seeds) currentRegions.add(seed.region);
        if (!currentRegions.equals(lastPreviewRegions))
            throw new IllegalArgumentException(
                "Lobule anchors changed after Preview.\n" +
                "Regenerate Preview before acceptance.");
        if (lastPreviewCoverage < 98.0)
            throw new IllegalArgumentException(String.format(Locale.US,
                "Acceptance stopped: only %.2f%% of the layer domain was assigned (minimum 98.00%%).\n" +
                "Fix guides/anchors and regenerate Preview.", lastPreviewCoverage));

        List<String> invalidCandidates = new ArrayList<String>();
        for (Integer boxed : previewIndexes) {
            Roi previewRoi = rois[boxed.intValue()];
            String previewName = safeName(previewRoi);
            String canonical = previewName.substring(AutoSplitEngine.AUTO_PREFIX.length());
            if (existing.contains(canonical) || !future.add(canonical)) collisions.add(canonical);
            String problem = candidateTypeProblem(canonical, previewRoi);
            if (problem != null) invalidCandidates.add(canonical + ": " + problem);
        }
        if (!invalidCandidates.isEmpty())
            throw new IllegalArgumentException(
                "Acceptance stopped because candidate names/types are invalid:\n- " +
                join(invalidCandidates, "\n- "));
        if (!collisions.isEmpty())
            throw new IllegalArgumentException(
                "Acceptance stopped to prevent overwrite/duplicate ROI names:\n- " + join(collisions, "\n- ") +
                "\n\nDelete or rename the conflicting final ROI, then review again.");

        List<String> missingCandidates = new ArrayList<String>();
        for (String region : seedRegions) {
            if (lastPreviewPL && !containsCandidatePrefix(future, "PL_" + region + "_part"))
                missingCandidates.add("PL_" + region);
            if (lastPreviewML && !containsCandidatePrefix(future, "ML_" + region + "_part"))
                missingCandidates.add("ML_" + region);
            if (lastPreviewGL && !containsCandidatePrefix(future, "GL_" + region + "_part"))
                missingCandidates.add("GL_" + region);
        }
        if (!missingCandidates.isEmpty())
            throw new IllegalArgumentException(
                "Acceptance stopped because expected candidates are missing:\n- " +
                join(missingCandidates, "\n- ") +
                "\n\nMove the relevant seed or add/extend a divider, then regenerate Preview.");

        String[] requiredTotals = {"Cb", "PL_TOTAL", "ML_TOTAL", "GL_TOTAL", "WM", "GL+WM"};
        List<String> missingTotals = new ArrayList<String>();
        for (String required : requiredTotals) {
            if (!existing.contains(required) && !future.contains(required)) missingTotals.add(required);
        }
        if (!missingTotals.isEmpty())
            throw new IllegalArgumentException(
                "Acceptance stopped because required whole/total ROIs would still be missing:\n- " +
                join(missingTotals, "\n- "));

        GenericDialog dialog = new GenericDialog(TITLE + " - Accept Human-reviewed Candidates");
        dialog.addMessage("Candidates: " + previewIndexes.size() + "\n" +
            "Accept removes AUTO_REVIEW_ from all candidates in one batch and applies normal CRA ROI names/colors.\n" +
            "You do not need to rename individual ROIs.\n" +
            "This action does not overwrite an existing final ROI.");
        dialog.addCheckbox("I visually checked every candidate boundary", false);
        dialog.setOKLabel("Accept candidates");
        dialog.showDialog();
        if (dialog.wasCanceled()) return;
        boolean visuallyChecked = dialog.getNextBoolean();
        if (!visuallyChecked)
            throw new IllegalArgumentException("Acceptance canceled: visual review confirmation was not checked.");

        String backupPath = savePreAcceptBackup(manager, image);

        for (Integer boxed : previewIndexes) {
            int index = boxed.intValue();
            Roi roi = manager.getRoi(index);
            String canonical = safeName(roi).substring(AutoSplitEngine.AUTO_PREFIX.length());
            Roi accepted = (Roi) roi.clone();
            accepted.setName(canonical);
            applyCanonicalColor(accepted, canonical);
            manager.setRoi(accepted, index);
            manager.rename(index, canonical);
        }
        deleteHelpers(false);
        image.deleteRoi();
        showAllWithLabels(manager, image);
        resetPreviewManifest();
        IJ.showMessage(TITLE + " - Accepted",
            previewIndexes.size() + " candidate ROI(s) accepted.\n" +
            "SEED_/DIVIDER_ helpers were removed from the working ROI Manager.\n" +
            "Pre-accept rollback ZIP: " + backupPath + "\n" +
            "FISSURE_ ROIs and original whole/total ROIs were preserved.\n\n" +
            "Run the existing CRA Structural QC / Final QC before export.");
    }

    private static int clearPreview(boolean notify) {
        RoiManager manager = RoiManager.getRoiManager();
        int removed = deleteByPrefixes(manager, new String[]{AutoSplitEngine.AUTO_PREFIX});
        resetPreviewManifest();
        ImagePlus image = currentImageOrNull();
        if (image != null) showAllWithLabels(manager, image);
        if (notify) IJ.showMessage(TITLE, removed + " preview candidate(s) removed. Original/final ROIs were preserved.");
        return removed;
    }

    private static String candidateTypeProblem(String canonical, Roi roi) {
        if (canonical.matches("PL_[A-Za-z0-9_.+\\-]+_part[1-9][0-9]*")) {
            if (!roi.isLine() || roi.isArea() || roi.getLength() <= 0.0)
                return "local PL must be a nonzero line/polyline ROI";
            return null;
        }
        if (canonical.matches("(ML|GL)_[A-Za-z0-9_.+\\-]+_part[1-9][0-9]*") ||
                "ML_TOTAL".equals(canonical) || "GL_TOTAL".equals(canonical)) {
            Rectangle bounds = roi.getBounds();
            if (!roi.isArea() || bounds.width <= 0 || bounds.height <= 0)
                return "ML/GL must be a nonzero area ROI";
            return null;
        }
        return "unrecognized AUTO_REVIEW name";
    }

    private static String savePreAcceptBackup(RoiManager manager, ImagePlus image) {
        String directory = IJ.getDirectory("image");
        if (directory == null || directory.length() == 0)
            throw new IllegalArgumentException(
                "Cannot determine the montage folder, so a rollback ZIP cannot be saved.\n" +
                "Save/open the montage from a folder and try again.");
        File folder = new File(directory);
        String sample = folder.getName();
        if (sample == null || sample.length() == 0) sample = "sample";
        String stamp = new SimpleDateFormat("yyyyMMdd_HHmmss", Locale.US).format(new Date());
        String path = new File(folder,
            sample + "_CRA_AUTOSPLIT_PREACCEPT_" + stamp + ".zip").getAbsolutePath();
        if (!manager.runCommand("Save", path) || !new File(path).isFile())
            throw new IllegalArgumentException("Failed to save the pre-accept rollback ZIP:\n" + path);
        return path;
    }

    private static void resetPreviewManifest() {
        lastPreviewNames.clear();
        lastPreviewRegions.clear();
        lastPreviewCoverage = -1.0;
        lastPreviewPL = false;
        lastPreviewML = false;
        lastPreviewGL = false;
    }

    private static int deleteHelpers(boolean notify) {
        RoiManager manager = RoiManager.getRoiManager();
        int removed = deleteByPrefixes(manager, new String[]{SEED_PREFIX, DIVIDER_PREFIX});
        ImagePlus image = currentImageOrNull();
        if (image != null) showAllWithLabels(manager, image);
        if (notify) IJ.showMessage(TITLE, removed + " SEED_/DIVIDER_ helper ROI(s) removed. FISSURE_ ROIs were preserved.");
        return removed;
    }

    private static Validation inspect(Roi[] rois) {
        Validation out = new Validation();
        Map<String, Integer> nameCounts = new LinkedHashMap<String, Integer>();
        for (Roi roi : rois) {
            String name = safeName(roi);
            if (name.isEmpty()) {
                out.errors.add("An unnamed ROI exists in ROI Manager.");
                continue;
            }
            Integer old = nameCounts.get(name);
            nameCounts.put(name, Integer.valueOf(old == null ? 1 : old.intValue() + 1));
            if (!out.namedRois.containsKey(name)) out.namedRois.put(name, (Roi) roi.clone());
            if (name.startsWith("FISSURE_")) out.fissureCount++;
            if (name.startsWith(DIVIDER_PREFIX)) out.dividerCount++;
        }
        for (Map.Entry<String, Integer> entry : nameCounts.entrySet()) {
            if (entry.getValue().intValue() > 1)
                out.errors.add("Duplicate ROI name: " + entry.getKey() + " (" + entry.getValue() + ")");
        }

        requireForValidation(out, "Cb", true);
        requireForValidation(out, "WM", true);
        requireForValidation(out, "GL+WM", true);
        requireForValidation(out, "PL_TOTAL", false);
        Roi plTotal = out.namedRois.get("PL_TOTAL");
        if (plTotal != null) {
            Rectangle bounds = plTotal.getBounds();
            if ((!plTotal.isArea() && !plTotal.isLine()) || bounds.width <= 0 || bounds.height <= 0)
                out.errors.add("PL_TOTAL must be a nonzero area band or line ROI, but is " +
                    plTotal.getTypeAsString() + ".");
        }
        if (!out.namedRois.containsKey("ML_TOTAL"))
            out.warnings.add("ML_TOTAL is absent; Preview can derive AUTO_REVIEW_ML_TOTAL.");
        if (!out.namedRois.containsKey("GL_TOTAL"))
            out.warnings.add("GL_TOTAL is absent; Preview can derive AUTO_REVIEW_GL_TOTAL.");

        Roi cb = out.namedRois.get("Cb");
        Roi wm = out.namedRois.get("WM");
        int explicitSeedCount = countByPrefix(rois, SEED_PREFIX);

        /* Existing named local PL ROIs are the least burdensome and most
         * traceable semantic source.  They may be legacy line ROIs OR traced
         * area ROIs; both are valid as name/position anchors and are preserved. */
        deriveSeedsFromLocalPL(rois, out);
        if (out.seeds.size() >= 2) {
            if (explicitSeedCount > 0)
                out.warnings.add(explicitSeedCount +
                    " SEED_ helper(s) were ignored because named local PL anchors are available.");
        } else {
            out.seeds.clear();
            out.derivedFromLocalPL = false;
            out.anchorSource = "";
            readExplicitSeeds(rois, out, cb, wm);
        }
        if (out.seeds.size() < 2) {
            if (explicitSeedCount == 0)
                out.errors.add("At least two lobule anchors are required. Add one-point SEED_ ROIs, " +
                    "or keep at least two verified local PL_<lobule>_part... lines in ROI Manager.");
            else
                out.errors.add("At least two valid single-point SEED_ ROIs are required.");
        }
        if (out.anchorSource.length() == 0)
            out.anchorSource = explicitSeedCount > 0 ? "Explicit one-point SEED_ ROIs" : "None";
        for (int i = 0; i < out.seeds.size(); i++) {
            AutoSplitEngine.Seed a = out.seeds.get(i);
            for (int j = i + 1; j < out.seeds.size(); j++) {
                AutoSplitEngine.Seed b = out.seeds.get(j);
                double dx = a.x - b.x;
                double dy = a.y - b.y;
                if (dx * dx + dy * dy < 9.0)
                    out.errors.add("Seeds " + a.region + " and " + b.region + " overlap.");
            }
        }
        if (out.fissureCount + out.dividerCount == 0)
            out.warnings.add("No FISSURE_ or DIVIDER_ guide exists; anatomical layer components may remain connected.");
        return out;
    }

    /** Reuse already verified legacy local PL lines or traced-area ROIs as
     * semantic lobule anchors.
     * This removes all SEED_ work in the recommended ML/GL-only workflow and
     * preserves the user's existing PL measurements unchanged. */
    private static void deriveSeedsFromLocalPL(Roi[] rois, Validation out) {
        Map<String, Roi> longestByRegion = new LinkedHashMap<String, Roi>();
        Map<String, Double> scoreByRegion = new LinkedHashMap<String, Double>();
        for (Roi roi : rois) {
            String name = safeName(roi);
            if (name.startsWith("REVIEW_")) name = name.substring("REVIEW_".length());
            if (name.startsWith(AutoSplitEngine.AUTO_PREFIX)) continue;
            if (!name.startsWith("PL_") || "PL_TOTAL".equals(name) ||
                    (!roi.isLine() && !roi.isArea())) continue;
            int part = name.lastIndexOf("_part");
            if (part <= 3) continue;
            String region = name.substring(3, part).trim();
            if (region.length() == 0) continue;
            double score = roi.isArea() ? maskPixelCount(roi) : roi.getLength();
            Double old = scoreByRegion.get(region);
            if (old == null || score > old.doubleValue()) {
                longestByRegion.put(region, roi);
                scoreByRegion.put(region, Double.valueOf(score));
            }
        }
        int lineCount = 0;
        int areaCount = 0;
        for (Map.Entry<String, Roi> entry : longestByRegion.entrySet()) {
            Roi source = entry.getValue();
            double[] point = representativePoint(source);
            if (point == null) continue;
            out.seeds.add(new AutoSplitEngine.Seed(entry.getKey(), point[0], point[1]));
            if (source.isArea()) areaCount++; else lineCount++;
        }
        if (!out.seeds.isEmpty()) {
            out.derivedFromLocalPL = true;
            out.anchorSource = "Existing local PL names/positions (line=" + lineCount +
                ", traced area=" + areaCount + "; no SEED_ required)";
            out.warnings.add("Local PL ROIs are used only as semantic anchors and are preserved unchanged.");
        }
    }

    private static void readExplicitSeeds(Roi[] rois, Validation out, Roi cb, Roi wm) {
        Set<String> seedRegions = new HashSet<String>();
        for (Roi roi : rois) {
            String name = safeName(roi);
            if (!name.startsWith(SEED_PREFIX)) continue;
            String region = name.substring(SEED_PREFIX.length()).trim();
            if (region.isEmpty()) {
                out.errors.add(name + " has no region name.");
                continue;
            }
            if (!(roi instanceof PointRoi)) {
                out.errors.add(name + " is not a Point ROI.");
                continue;
            }
            FloatPolygon point = roi.getFloatPolygon();
            int count = point == null ? 0 : point.npoints;
            if (count != 1) {
                out.errors.add(name + " contains " + count + " points; exactly one is required. " +
                    "Re-click the intended lobule and use Add Lobule Seed.");
                continue;
            }
            if (!seedRegions.add(region)) {
                out.errors.add("Duplicate seed region: " + region);
                continue;
            }
            double x = point.xpoints[0];
            double y = point.ypoints[0];
            if (cb != null && cb.isArea() && !cb.containsPoint(x, y))
                out.errors.add(name + " lies outside Cb.");
            if (wm != null && wm.isArea() && wm.containsPoint(x, y))
                out.warnings.add(name + " lies inside WM; move it into the lobule cortical ribbon for a more stable split.");
            out.seeds.add(new AutoSplitEngine.Seed(region, x, y));
        }
    }

    private static int countByPrefix(Roi[] rois, String prefix) {
        int count = 0;
        for (Roi roi : rois) if (safeName(roi).startsWith(prefix)) count++;
        return count;
    }

    /** Finds an interior representative point without changing the ROI. */
    private static double[] representativePoint(Roi roi) {
        if (roi.isLine()) {
            FloatPolygon points = roi.getFloatPolygon();
            if (points == null || points.npoints < 1) return null;
            int middle = points.npoints / 2;
            return new double[]{points.xpoints[middle], points.ypoints[middle]};
        }
        Rectangle bounds = roi.getBounds();
        ImageProcessor mask = roi.getMask();
        if (mask == null || mask.getWidth() <= 0 || mask.getHeight() <= 0) {
            double x = bounds.getCenterX();
            double y = bounds.getCenterY();
            return roi.containsPoint(x, y) ? new double[]{x, y} : null;
        }
        long count = 0;
        double sumX = 0.0;
        double sumY = 0.0;
        for (int y = 0; y < mask.getHeight(); y++) {
            for (int x = 0; x < mask.getWidth(); x++) {
                if (mask.get(x, y) == 0) continue;
                count++;
                sumX += x;
                sumY += y;
            }
        }
        if (count == 0) return null;
        int cx = (int) Math.round(sumX / count);
        int cy = (int) Math.round(sumY / count);
        int bestX = -1;
        int bestY = -1;
        long bestDistance = Long.MAX_VALUE;
        for (int y = 0; y < mask.getHeight(); y++) {
            for (int x = 0; x < mask.getWidth(); x++) {
                if (mask.get(x, y) == 0) continue;
                long dx = x - cx;
                long dy = y - cy;
                long distance = dx * dx + dy * dy;
                if (distance < bestDistance) {
                    bestDistance = distance;
                    bestX = x;
                    bestY = y;
                }
            }
        }
        if (bestX < 0) return null;
        return new double[]{bounds.x + bestX + 0.5, bounds.y + bestY + 0.5};
    }

    private static double maskPixelCount(Roi roi) {
        ImageProcessor mask = roi.getMask();
        if (mask == null) {
            Rectangle bounds = roi.getBounds();
            return Math.max(1.0, (double) bounds.width * bounds.height);
        }
        long count = 0;
        for (int y = 0; y < mask.getHeight(); y++)
            for (int x = 0; x < mask.getWidth(); x++)
                if (mask.get(x, y) != 0) count++;
        return Math.max(1.0, count);
    }

    private static void requireForValidation(Validation validation, String name, boolean mustBeArea) {
        Roi roi = validation.namedRois.get(name);
        if (roi == null) {
            validation.errors.add("Required ROI is missing: " + name);
        } else if (mustBeArea && !roi.isArea()) {
            validation.errors.add(name + " must be an area ROI, but is " + roi.getTypeAsString() + ".");
        }
    }

    private static String validationFailureMessage(Validation validation) {
        StringBuilder message = new StringBuilder("Preview cannot start. Fix these items first:\n\n");
        for (String error : validation.errors) message.append("- ").append(error).append('\n');
        if (!validation.warnings.isEmpty()) {
            message.append("\nWarnings:\n");
            for (String warning : validation.warnings) message.append("- ").append(warning).append('\n');
        }
        return message.toString();
    }

    private static Map<String, Roi> cloneNamedRois(Map<String, Roi> source) {
        Map<String, Roi> copy = new LinkedHashMap<String, Roi>();
        for (Map.Entry<String, Roi> entry : source.entrySet())
            copy.put(entry.getKey(), (Roi) entry.getValue().clone());
        return copy;
    }

    private static void applyCanonicalColor(Roi roi, String name) {
        roi.setFillColor(null);
        if (name.startsWith("PL_") || "PL_TOTAL".equals(name)) {
            roi.setStrokeColor(Color.RED);
        } else if (name.startsWith("ML_") || "ML_TOTAL".equals(name)) {
            roi.setStrokeColor(new Color(0, 180, 0));
        } else if (name.startsWith("GL_") || "GL_TOTAL".equals(name)) {
            roi.setStrokeColor(Color.YELLOW);
        } else if (name.startsWith("FISSURE_")) {
            roi.setStrokeColor(new Color(255, 153, 0));
        } else {
            roi.setStrokeColor(Color.WHITE);
        }
        roi.setStrokeWidth(2.0);
    }

    private static ImagePlus requireImage() {
        ImagePlus image = currentImageOrNull();
        if (image == null) throw new IllegalArgumentException("Open a montage image first.");
        return image;
    }

    private static ImagePlus currentImageOrNull() {
        try {
            return IJ.getImage();
        } catch (RuntimeException noImage) {
            return null;
        }
    }

    private static void showAllWithLabels(RoiManager manager, ImagePlus image) {
        manager.runCommand(image, "Show All with labels");
        image.updateAndDraw();
    }

    private static int deleteByPrefixes(RoiManager manager, String[] prefixes) {
        int removed = 0;
        for (int i = manager.getCount() - 1; i >= 0; i--) {
            String name = manager.getName(i);
            for (String prefix : prefixes) {
                if (name != null && name.startsWith(prefix)) {
                    deleteIndex(manager, i);
                    removed++;
                    break;
                }
            }
        }
        return removed;
    }

    private static void deleteExactName(RoiManager manager, String target) {
        for (int i = manager.getCount() - 1; i >= 0; i--) {
            if (target.equals(manager.getName(i))) deleteIndex(manager, i);
        }
    }

    private static void deleteIndex(RoiManager manager, int index) {
        manager.setSelectedIndexes(new int[]{index});
        manager.runCommand("Delete");
    }

    private static String safeName(Roi roi) {
        return roi == null || roi.getName() == null ? "" : roi.getName().trim();
    }

    private static String sanitizeRegion(String source) {
        String text = source == null ? "" : source.trim();
        if (text.isEmpty()) throw new IllegalArgumentException("Lobule/divider name cannot be empty.");
        text = text.replace("4/5Cb", "4-5Cb");
        text = text.replace('/', '-').replace(' ', '_');
        text = text.replaceAll("[^A-Za-z0-9_.+-]", "-");
        while (text.contains("--")) text = text.replace("--", "-");
        if (text.isEmpty()) throw new IllegalArgumentException("The name contains no usable character.");
        return text;
    }

    private static void appendMessages(StringBuilder target, String label, List<String> messages) {
        if (messages.isEmpty()) return;
        target.append('\n').append(label).append("S:\n");
        for (String message : messages) target.append("- ").append(message).append('\n');
    }

    private static String join(List<String> values, String separator) {
        StringBuilder result = new StringBuilder();
        for (int i = 0; i < values.size(); i++) {
            if (i > 0) result.append(separator);
            result.append(values.get(i));
        }
        return result.toString();
    }

    private static boolean containsCandidatePrefix(Set<String> canonicalNames, String prefix) {
        for (String name : canonicalNames) if (name.startsWith(prefix)) return true;
        return false;
    }

    private static void showGuide() {
        IJ.showMessage(TITLE + " - Safe workflow",
            "1. Keep verified Cb, PL_TOTAL, WM, GL+WM and FISSURE_ ROIs.\n" +
            "2A. To auto-generate PL too: Multi-point tool -> click each lobule once in order -> Add Seed Set.\n" +
            "    With SEED_ anchors and a line PL_TOTAL, enable \"Create PL candidates\" in Preview\n" +
            "    to split PL_TOTAL's own traced vertices per lobule (exact length, no manual per-lobule PL).\n" +
            "2B. If you already have verified local PL_<lobule>_part... ROIs and want to keep them\n" +
            "    unchanged instead: skip SEED_, keep \"Create PL candidates\" OFF, and only ML/GL\n" +
            "    are generated (local PL is reused only as the lobule anchor).\n" +
            "3. Run Smart Start. It validates inputs and opens Preview only when ready.\n" +
            "4. AUTO_REVIEW_ ROIs are proposals only (PL red, ML green, GL yellow).\n" +
            "5. If a split is wrong, draw a line across that boundary and Add Divider, then Preview again.\n" +
            "6. Inspect every PL/ML/GL boundary visually.\n" +
            "7. Accept only after visual review. Existing final ROI names are never overwritten.\n" +
            "8. Accept removes AUTO_REVIEW_ in one batch; no manual renaming is needed.\n" +
            "   Temporary SEED_/DIVIDER_ helpers are deleted by default.\n" +
            "9. Keep using the existing CRA Structural QC, Final QC and export functions.\n\n" +
            "Important: legacy SEED_ ROIs that contain accumulated points are invalid and must be re-added.");
    }

    private static final class Validation {
        final Map<String, Roi> namedRois = new LinkedHashMap<String, Roi>();
        final List<AutoSplitEngine.Seed> seeds = new ArrayList<AutoSplitEngine.Seed>();
        final List<String> errors = new ArrayList<String>();
        final List<String> warnings = new ArrayList<String>();
        int fissureCount;
        int dividerCount;
        String anchorSource = "";
        boolean derivedFromLocalPL;
    }
}
