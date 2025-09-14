네, Range Shifter를 사용할 때마다 실시간으로 계산하는 방식의 비효율성에 대한 지적은 매우 정확합니다. 환자 개개인의 계산 속도를 극대화하려는 MOQUI의 목적에 맞게, **Range Shifter의 효과를 미리 계산(pre-calculation)하여 테이블 형태로 만들어두고, 시뮬레이션 시에는 이 테이블을 참조**하는 방식이 훨씬 효율적입니다.

TOPAS를 이용해 다양한 Range Shifter 위치에 따른 위상 공간(Phase Space) 정보를 생성하고, 이 데이터를 MOQUI에서 사용할 수 있는 테이블 형태로 통합하는 전체적인 워크플로우와 코드 수정 방향을 단계별로 상세히 안내해 드리겠습니다.

-----

### \#\# 📝 1단계: TOPAS를 이용한 위상 공간 데이터 생성

첫 번째 단계는 TOPAS를 사용하여 'Range Shifter 바로 다음' 지점에서의 입자 정보를 생성하는 것입니다. 이 작업을 **모든 빔 에너지**와 **사용 가능한 모든 Range Shifter 위치/두께 조합**에 대해 반복해야 합니다.

1.  **TOPAS 시뮬레이션 환경 구축**:

      * 실제 치료 기기의 노즐을 Geant4/TOPAS 지오메트리로 최대한 정확하게 모델링합니다. (빔 소스부터 Range Shifter가 위치할 수 있는 모든 지점까지)
      * `TsPhaseSpace` Scorer를 사용하여 위상 공간 정보를 저장할 가상의 평면(scoring plane)을 생성합니다. 이 평면은 **각 Range Shifter 바로 뒷단**에 위치해야 합니다.

2.  **TOPAS 파라미터 파일 (`.topas`) 설정**:

      * Scorer 설정을 통해 필요한 모든 정보를 기록하도록 지정합니다. 최소한 다음 7가지 정보는 필수적입니다.
          * X, Y 위치 (Position)
          * X', Y' 방향 (Direction Cosines, 보통 `Px/Pz`, `Py/Pz`로 계산)
          * 운동에너지 (Kinetic Energy)
          * 입자 종류 (Particle Type, 양성자인지 확인)
          * 가중치 (Weight, 보통 1)
      * 아래는 TOPAS 파라미터 파일의 예시입니다.

    <!-- end list -->

    ```topas
    # --- Scorer Definition ---
    s:Sc/MyPhaseSpace/Type = "PhaseSpace"
    s:Sc/MyPhaseSpace/Surface = "RS_ExitPlane"  # Range Shifter 바로 뒤에 위치한 가상 평면
    s:Sc/MyPhaseSpace/OutputFile = "./phsp/E100_RS_pos1" # 에너지와 RS 위치 조합별로 파일명 지정
    s:Sc/MyPhaseSpace/OutputFormat = "ASCII" # 또는 "Binary"
    b:Sc/MyPhaseSpace/IncludeBeamName = "False"

    # 기록할 변수들 정의
    sv:Sc/MyPhaseSpace/Report = 8 "PositionX" "PositionY" "DirectionX" "DirectionY" "Energy" "ParticleName" "Weight" "FlagIsNewHistory"
    ```

3.  **시뮬레이션 실행 및 데이터 생성**:

      * 모든 `에너지-RS 위치` 조합에 대해 TOPAS 시뮬레이션을 실행하여 위상 공간 파일(`.phsp`)들을 생성합니다. 이 파일들은 수백만 개 입자의 개별 정보가 담긴 거대한 데이터셋이 됩니다.

-----

### \#\# 📊 2단계: 위상 공간 파일 분석 및 파라미터화

수십, 수백 개의 거대한 위상 공간 파일을 MOQUI에서 직접 읽는 것은 비효율적입니다. MOQUI의 강점은 통계적 모델을 사용하는 것이므로, 이 위상 공간 파일들을 분석하여 MOQUI가 사용하는 **트위스 파라미터(Twiss Parameters)와 에너지 분포로 다시 파라미터화**해야 합니다.

1.  **분석 스크립트 작성**: Python (uproot, pandas 라이브러리) 또는 C++/ROOT를 사용하여 각 위상 공간 파일을 읽고 분석하는 스크립트를 작성합니다.
2.  **통계적 파라미터 계산**: 각 파일에 대해 다음 값들을 계산합니다.
      * **에너지**: 평균 에너지 (`<E>`), 에너지 표준편차 (`σE`)
      * **X 평면 트위스 파라미터**: `<x>`, `<x'>`, `<x²>`, `<x'²>`, `<xx'>` 값을 계산하여 αx, βx, εx (emittance)를 구합니다.
      * **Y 평면 트위스 파라미터**: `<y>`, `<y'>`, `<y²>`, `<y'²>`, `<yy'>` 값을 계산하여 αy, βy, εy를 구합니다.
3.  **데이터 테이블 생성**: 분석된 결과를 다음과 같은 형태의 테이블(예: CSV 파일)로 정리합니다.

| Initial\_Energy (MeV) | RS\_ID | Mean\_Energy (MeV) | Sigma\_Energy (MeV) | Alpha\_X | Beta\_X (mm/mrad) | Emit\_X (mm*mrad) | Alpha\_Y | Beta\_Y (mm/mrad) | Emit\_Y (mm*mrad) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 100.0 | 1 | 95.2 | 0.8 | 0.12 | 2.5 | 3.1 | -0.05 | 2.8 | 2.9 |
| 100.0 | 2 | 90.5 | 1.2 | 0.15 | 3.1 | 3.5 | -0.08 | 3.3 | 3.2 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| 230.0 | 3 | 221.0 | 0.6 | 0.08 | 1.8 | 2.5 | -0.03 | 1.9 | 2.4 |

-----

### \#\# 💻 3단계: MOQUI 코드 수정 및 테이블 통합

마지막으로, 생성된 데이터 테이블을 MOQUI 코드에 통합하고, Range Shifter가 사용될 때 이 테이블을 참조하도록 로직을 수정합니다.

1.  **데이터 테이블 내장**:

      * 2단계에서 생성한 테이블을 `treatment_machines/` 폴더 안의 기기 전용 헤더 파일(예: **`mqi_treatment_machine_smc_gtr2.hpp`**)에 **`static const` C++ 배열** 형태로 변환하여 직접 삽입합니다. 이는 외부 파일 I/O 없이 가장 빠르게 데이터에 접근하는 방법입니다.

    <!-- end list -->

    ```cpp
    // In mqi_treatment_machine_smc_gtr2.hpp
    namespace mqi {
        // 새로운 데이터 구조체 정의
        struct rs_phsp_parameters {
            float initial_energy;
            int   rs_id;
            float mean_energy;
            // ... all other twiss parameters
        };

        // TOPAS로 계산된 데이터를 static 배열로 선언
        static const rs_phsp_parameters gtr2_rs_phsp_table[] = {
            { 100.0f, 1, 95.2f, 0.8f, ... },
            { 100.0f, 2, 90.5f, 1.2f, ... },
            // ... all data
        };
    }
    ```

2.  **빔 소스 생성 로직 수정 (`mqi_treatment_machine_ion.hpp`, `mqi_beamsource.hpp`):**

      * **`mqi::treatment_machine_ion<T>::create_beamsource(...)`** 함수를 수정합니다.
      * 이 함수는 DICOM 데이터를 파싱하여 `beamsource` 객체를 생성합니다. 이 때, 치료 계획에 Range Shifter가 사용되었는지 (`RangeShifterID` 태그 확인)를 검사하는 로직을 추가합니다.
      * **If (Range Shifter is used):**
          * `gtr2_rs_phsp_table`에서 현재 스팟의 초기 에너지와 RS ID에 해당하는 파라미터 세트를 검색하는 함수를 호출합니다.
          * 검색된 **'RS 통과 후'** 파라미터(새로운 평균 에너지, 트위스 파라미터 등)를 `mqi::beamlet` 생성자에 전달하여 입자를 샘플링하도록 합니다.
      * **Else (No Range Shifter):**
          * 기존과 동일하게 '공기 중(in-air)' 위상 공간 파라미터를 사용하여 입자를 샘플링합니다.

3.  **빔라인(Beamline) 수정 (`mqi_treatment_machine_ion.hpp`):**

      * **`mqi::treatment_machine_ion<T>::create_beamline(...)`** 함수를 수정합니다.
      * 여기서는 반대로, 위상 공간 테이블을 사용하기로 결정했다면 **`beamline`에 `mqi::rangeshifter` 객체를 추가하는 로직을 비활성화**해야 합니다. 그렇지 않으면 Range Shifter의 효과가 중복으로 계산됩니다.

이러한 **'선계산(pre-calculation) 및 테이블 참조(look-up table)'** 방식은 초기 설정에 품이 많이 들지만, 일단 구축되고 나면 MOQUI의 가장 큰 장점인 **실시간 계산 속도를 희생하지 않으면서도 Range Shifter의 물리적 효과를 정확하게 반영**할 수 있는 가장 이상적인 접근법입니다.