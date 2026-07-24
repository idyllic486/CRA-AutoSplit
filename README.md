# CRA AutoSplit

Fiji/ImageJ1 plugin that auto-generates per-lobule **PL / ML / GL** ROI
candidates for cerebellum montage annotation, as a review-and-accept
assistant on top of the existing **CRA Legacy** toolbar/macro workflow
(see [FIJI-Macro](https://github.com/idyllic486/FIJI-Macro)).

> 기존 **CRA Legacy** 워크플로를 대체하지 않고, PL/ML/GL local ROI를 손으로
> 그리는 시간을 줄여주는 보조 플러그인입니다. 자동 생성 결과는 최종 판정이
> 아니라 **검토용 후보**(`AUTO_REVIEW_` 접두사)이며, 연구자가 시각적으로
> 확인하고 승인(Accept)해야 최종 ROI 이름으로 확정됩니다.

## 요구사항

- Fiji (ImageJ 1.54p 이상), macOS/Windows/Linux
- Java 8 이상 (Fiji 번들 JDK로 충분, 별도 설치 불필요)

## 설치

1. Fiji를 완전히 종료합니다.
2. 이 저장소의 [`CRA_AutoSplit_A003_alpha4.jar`](CRA_AutoSplit_A003_alpha4.jar)를
   `Fiji.app/plugins/`에 복사합니다. (이전 alpha 버전의 jar가 있다면 함께 있어도
   되지만, 같은 Java 클래스를 담은 jar 두 개를 동시에 두면 클래스로더 충돌이
   날 수 있으므로 옛 버전은 지우거나 이 파일로 덮어쓰는 것을 권장합니다.)
3. `SHA256.txt`로 다운로드한 jar의 무결성을 확인할 수 있습니다:
   ```bash
   shasum -a 256 -c SHA256.txt
   ```
4. Fiji를 다시 실행합니다. `Plugins > CRA AutoSplit A003 alpha4` 메뉴가 보이면 설치 완료입니다.

## 빠른 시작

### A. 전체 자동화 (PL 포함) — `SEED_` 앵커 방식

1. 기존 CRA로 `Cb`, `PL_TOTAL`(선), `WM`, `GL+WM`, `FISSURE_*`를 준비합니다.
2. Multi-point 도구로 소엽마다 한 번씩 해부학적 순서대로 클릭 → **"1A. Add Seed Set"**
   실행 → 대화상자에 클릭 순서대로 소엽 이름을 쉼표로 입력 (예: `Sim,Crus1,Crus2,PM,PFl,Fl`).
3. **"3. Validate Inputs"** 실행 → ERROR가 없는지 확인.
4. **"4. Generate Exclusive-layer Preview"** 실행 → 설정창에서 **PL/ML/GL 후보 생성**을
   모두 켭니다. (`PL_TOTAL`이 선이면 "Create PL candidates (direct split of
   PL_TOTAL's traced vertices)" 옵션이 원본 좌표를 그대로 사용해 소엽별로 잘라줍니다.)
5. PL(빨강)·ML(초록)·GL(노랑) 경계를 전부 시각적으로 확인합니다.
6. 잘못 나뉜 경계가 있으면 그 경계를 가로지르는 선을 그리고 **"2. Add Divider Guide"**
   추가 후 Preview를 다시 실행합니다.
7. 모든 후보가 맞으면 **"5. Accept and Bulk-finalize Names"** 실행.
8. 기존 CRA Legacy의 Structural QC / Final QC / 결과 파일 생성 기능을 이어서 사용합니다.

### B. ML/GL만 자동화 — 기존 local PL 보존 방식

이미 소엽별로 검증된 `PL_<lobule>_partN` ROI(선 또는 traced area)가 있다면, `SEED_`
없이 그 이름/위치를 앵커로 그대로 재사용할 수 있습니다. 이 경우 Preview에서
PL 생성 옵션은 반드시 꺼둔 채(자동으로 꺼짐) ML/GL만 생성됩니다 — 기존에 확정한
PL 트레이스를 손대지 않습니다.

## 주요 기능

- PL/ML/GL 마스크를 상호 배타적으로 생성 (겹침 0 픽셀을 내부적으로 검증)
- **inner-ML sliver 자동 제거**: PL 안쪽(GL+WM과 PL 사이)에 갇히는 얇은 ML
  조각을 pial-surface flood-fill 방식으로 제거해, 소엽당 ML이 여러 조각으로
  쪼개지지 않고 하나의 깔끔한 바깥쪽 band로 나옵니다.
- **PL 직접 분할**: `PL_TOTAL`을 그릴 때 찍은 원본 좌표를 그대로 소엽별로
  잘라서 사용 — 래스터화/스켈레톤화로 인한 길이 손실이 없습니다.
- FISSURE_/DIVIDER_ 가이드 기반 경계 제어, 최소 컴포넌트 면적 필터링
- `AUTO_REVIEW_` 임시 접두사 + 승인 전 시각 검토 강제, 승인 시 롤백 ZIP 자동 백업
- 기존 CRA Legacy ROI 명명 규칙과 100% 호환

## 알려진 제한

- Fissure가 실제 layer를 충분히 가르지 못하면 일부 lobule 후보가 생성되지
  않을 수 있습니다 (Divider 추가 필요).
- PL 직접 분할은 `SEED_` 앵커 + line `PL_TOTAL` 조합에서만 적용됩니다. 기존
  local PL을 anchor로 재사용하는 워크플로(위 B)에서는 검증된 PL을 보존하기
  위해 계속 꺼져 있습니다.
- `PL_TOTAL`이 area/band ROI인 경우 PL 후보 생성은 여전히 실험적(스켈레톤
  기반) 방식이며 기본값이 꺼짐(off)입니다.
- 개발 검증은 제한된 수의 montage 샘플을 기준으로 했습니다 (`VALIDATION/` 참고).
  독립적인 대규모 생물학적 정확도 검증은 아닙니다.

## 검증

[`VALIDATION/`](VALIDATION/) 폴더에 입력 호환성 테스트와 수동 완료본 대비 비교
결과(Dice 계수, PL 길이 비율 등)가 있습니다.

## 소스에서 빌드하기

Maven/Ant 없이 Fiji에 내장된 JDK로 바로 빌드할 수 있습니다:

```bash
FIJI=/path/to/Fiji.app   # 또는 Fiji 실행 파일이 있는 폴더
JAVAC="$FIJI/java/*/*/*/bin/javac"   # 실제 경로는 환경마다 다름 — Finder에서 Fiji/java 폴더 확인
IJ_JAR="$FIJI/jars/ij-*.jar"

mkdir -p build/classes
javac -encoding UTF-8 --release 8 -cp $IJ_JAR -d build/classes SOURCE/*.java
cp SOURCE/plugins.config build/classes/
cd build/classes && jar cf ../CRA_AutoSplit_A003_alpha4.jar plugins.config cra
```

빌드된 jar를 `Fiji.app/plugins/`에 넣고 Fiji를 재시작하면 반영됩니다.

## 라이선스

[MIT](LICENSE)

## 변경 이력

[CHANGELOG.md](CHANGELOG.md)

## 관련 저장소

- [FIJI-Macro](https://github.com/idyllic486/FIJI-Macro) — CRA Legacy 8버튼
  툴바 + 단축키 매크로 (`.ijm`). 이 플러그인이 함께 동작을 전제로 하는
  ROI 명명 규칙과 워크플로가 정의된 곳입니다.
