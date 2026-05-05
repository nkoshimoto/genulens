# genulens status

## 最終目標

genulensのシミュレーション本体をオブジェクト化し、以下を実現する:

1. **Pythonから直接呼び出せる構造** — stdoutキャプチャ経由ではなく、EventSampler が Event/SimulationResult を直接返す
2. **CLI と C++ API の分離** — CLI は表示層、C++ API は typed result を返す層にする
3. **active_state マクロの廃止** — DensityModel / KinematicModel / StellarPopulation などの
   明示的なクラスに置き換え、RunContext へのグローバルポインタに依存しない設計にする
4. **pre_gapmoe との統一** — gapmoe_src/ の3ツールを genulens の共有モデルクラスで書き直し、
   同一の GalacticModel / EventSampler を使い回す

現時点の優先順位は 1 と 2。`pre_gapmoe` と packaging は後で扱う。

## 現時点の構造

### 実行フロー

```
genulens.cpp
  -> simulation/cli.cpp::run_cli()
      -> Initializer::create_context()
      -> Initializer::initialize_rng()
      -> Sampler::run_cli()
          -> Initializer::read_model_options()
          -> Initializer::finalize_spatial_model()
          -> Initializer::read_sampling_options()
          -> PopulationRuntime (IMF / LF テーブル)
          -> KinematicRuntimeTables (Shu DF テーブル)
          -> NsdMomentRuntime (NSD モーメントテーブル)
          -> ObservationConfig::from_cli()
          -> LineOfSightDensityGrid::build()
          -> MassFunction::init_from_population()
          -> EventSampler::run_cli()   <- MCループ本体
              -> CliEventReporter      <- MC行出力・summary出力
          -> cleanup
```

Python は `SimulationBackend::simulate()` から `Sampler::simulate()` を呼び、
`EventSampler` の event sink で `SimulationResult` を直接構築する。
以前の stdout キャプチャ + `VERBOSITY=3` パース経路は廃止済み。
Python 公開 `Config` の主要フィールド (`seed`, `l`, `b`, `n_simu`, `observed_tE`,
`observed_tE_error`, `model.imf`) は typed config として `Sampler` に渡す。
加えて `Config.observation`, `Config.source`, `Config.sampling`, `Config.runtime` を追加し、
観測制約、source/extinction、sampling 制御、距離グリッドの基本設定も typed にした。
`Config.model` も `imf`, `density`, `kinematics`, `nsd`, `bh_kick` に分割し、
主要な銀河モデル・運動学・NSD・BH kick 設定を Python からオブジェクトとして指定できる。
`raw_cli_args` はまだ typed 化していない legacy option 用の互換口として残す。

`Sampler` はまだ大きいが、`PreparedSimulation` を独立した simulation object として導入し、
population/kinematic tables/NSD moments/grid/mass function/event config/observation を
準備済み状態として束ねた。イベント実行は `run_prepared_events()` に分離し、
`EventSampler::run()` は CLI 非依存の実行口になった。
MC 行出力と summary 出力は `CliEventReporter` に分離済み。

### 目標フロー

```
Python / C++ API
  -> EventSimulator::simulate()
      -> SimulationBackend::simulate()
          -> Initializer + Sampler preparation
          -> EventSampler::run(..., on_event, custom_likelihood)
              -> SimulationResult.events に Event を追加

CLI
  -> run_cli()
      -> 同じ preparation
      -> EventSampler::run(..., on_event)
          -> CliEventReporter で表示
```

重要な境界:

- `EventSampler` はイベントを生成し、観測尤度と custom likelihood をサンプリングループ内で適用する。
- CLI の MC 出力・summary は `CliEventReporter` に閉じ込める。
- Python API は `VERBOSITY=3` の文字列を読まない。
- `SimulationResult` は `Event` の typed collection として ndarray 変換を担当する。

### 主要ファイル (simulation/)

| ファイル | 役割 |
|---|---|
| sampler.cpp | オーケストレーター。まだ大きく、準備処理の分割が残る |
| prepared_simulation.{hpp,cpp} | 準備済み runtime state と cleanup/event 実行境界 |
| los_density_grid.{hpp,cpp} | 視線密度グリッド構築・サンプリング・BH補正 |
| observation_config.{hpp,cpp} | 観測制約 + 重要度サンプリング範囲 |
| mass_function.{hpp,cpp} | レンズ質量サンプリング・残骸進化・等級変換 |
| event_sampler.{hpp,cpp} | MCループ・統計更新・EventSink通知 |
| event_reporter.{hpp,cpp} | CLI の Monte Carlo 行出力・summary出力 |
| velocity_distribution.{hpp,cpp} | 速度サンプリング薄ラッパー |
| initialize.{hpp,cpp} | モデルオプション読み込み・RunContext初期化 |
| run_context.hpp | シミュレーション可変状態の集約 |
| internal/runtime.hpp | active_state宣言・互換マクロ (過渡的) |

### 互換ブリッジ (active_state)

`extern RunContext *active_state` + マクロ群が全 runtime .cpp ファイルで生きている。
各メソッドの先頭で `active_state = &ctx` をセットすることで科学計算部分の式は変えずに動く。
これは次の廃止ターゲット。

## 次のステップ

1. Sampler の準備処理をさらに分割する — option/config 適用、runtime table 準備、grid 準備を別関数/クラスにする
2. Sampler 内に残る CLI header/input parameter printf を reporter/formatter 側へ移す
3. `raw_cli_args` 互換口を縮小する — Python API では typed option を優先し、CLI option 文字列を隠蔽する
4. Config の責務整理 — `GenulensConfig` が大きくなってきたため、公開 API と内部 runtime config の境界を明確にする
5. active_state マクロの廃止 — runtime 関数を明示的参照を取るクラスに変換
6. pre_gapmoe rewrite — gapmoe_src/ の書き直し (.gaplan/gapmoe_rewrite.md 参照)
