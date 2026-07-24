# CRA AutoSplit A003 alpha4

## alpha4 — inner-ML sliver fix & PL auto-split

- **Fixed:** the inner-ML sliver (a thin strip of tissue trapped between PL
  and the GL+WM boundary) was not actually being removed by either ML
  derivation path, despite the original alpha4 design intending to remove it.
  Root cause: `PL_TOTAL` is traced as a single continuous **open** curve
  across the whole cortical ribbon (it runs through every fissure rather than
  closing into a separate loop per lobule). An open arc never separates the
  plane, so the previous `fillHoles()` border-flood approach could always
  sneak around through the trace's two far-away open endpoints and reach
  every lobule's interior pocket — it silently filled nothing.
  Replaced with: seal PL_TOTAL's two true open endpoints to GL+WM, then keep
  only the candidate ML pixels reachable from the pial (outer Cb) surface
  without crossing PL (`pialSideMlMask`). Applied uniformly whether `ML_TOTAL`
  is supplied or derived.
- **Added:** direct per-lobule PL splitting (`splitPLDirect`). When
  `PL_TOTAL` is a line and lobule anchors come from `SEED_` points, PL is now
  cut into `PL_<lobule>_partN` candidates directly from `PL_TOTAL`'s own
  traced vertices — using the same barrier + nearest-seed ownership as
  ML/GL — instead of rasterizing and skeletonizing. This preserves the exact
  hand-traced PL length and removes the need to hand-draw local PL per
  lobule when using the `SEED_`-anchored workflow. Falls back to the existing
  skeleton-centreline method (still labeled EXPERIMENTAL) when `PL_TOTAL` is
  an area/band instead of a line.
- Removed dead/unused code left over from an earlier, never-wired-up
  grid-based attempt at PL splitting (`Options.gridScale`,
  `rasterizeSeparators`, `drawThickSegment`, `findNearestValid`,
  `nearbyLabel`, `labelForCoordinate`, `createLineCandidates`).

## alpha4 (original)

- Changed ML derivation from plain `Cb - (GL+WM) - PL` to a PL-envelope
  method that removes closed PL interiors and false inner-ML slivers.
- A supplied area `ML_TOTAL` now has precedence and is clipped to Cb while
  remaining exclusive from PL and GL+WM.
- ML topology and measurement masks now use the same molecular-side mask,
  preventing boundary mismatch fragments from becoming extra `part` ROIs.

## alpha3

## Added

- Smart Start: validation 통과 시에만 Preview로 진행
- Traced area local PL anchor 지원
- Line/area anchor 개수 표시
- local PL이 존재할 때 stale SEED 무시
- 입력검사 메시지 단순화

## Preserved

- mutually exclusive PL/ML/GL masks
- bulk removal of `AUTO_REVIEW_`
- pre-accept rollback ZIP
- helper cleanup
- candidate completeness and collision checks
- 98% minimum assigned-domain acceptance gate
- CRA Legacy post-QC/export workflow

## Not adopted

- 연결 layer를 강제로 Voronoi 분할하는 실험 변경은 50M-43에서 ML Dice가 낮아져 제외했습니다.
- 자동 후보 수를 늘리는 것보다 기존 경계 정확도를 우선했습니다.

## Known limitation

- 불완전한 fissure로 같은 connected component에 남아 있는 lobule은 Divider가 필요할 수 있습니다.
- PL 자동 분할은 `SEED_` 앵커 + line `PL_TOTAL` 조합에서만 새 direct-vertex 방식이 적용됩니다.
  기존 local PL을 anchor로 재사용하는 워크플로에서는 계속 꺼져 있습니다(검증된 PL을 보존하기 위함).
