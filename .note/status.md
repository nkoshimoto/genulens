# genulens status

## 最終目標

genulensのシミュレーション本体をオブジェクト化し、以下を実現する:

1. **Pythonから直接呼び出せる構造** — stdoutキャプチャ経由ではなく、EventSamplerが SimulationResult を直接返す
2. **active_state マクロの廃止** — DensityModel / KinematicModel / StellarPopulation などの
   明示的なクラスに置き換え、RunContext へのグローバルポインタに依存しない設計にする
3. **pre_gapmoe との統一** — gapmoe_src/ の3ツールを genulens の共有モデルクラスで書き直し、
   同一の GalacticModel / EventSampler を使い回す

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
          -> cleanup
```

Python は同じ Sampler::run_cli() を VERBOSITY=3 で叩き、stdout をキャプチャしてパースする
(SimulationBackend)。これが次の廃止対象。

### 主要ファイル (simulation/)

| ファイル | 役割 |
|---|---|
| sampler.cpp | 薄いオーケストレーター (~350行) |
| los_density_grid.{hpp,cpp} | 視線密度グリッド構築・サンプリング・BH補正 |
| observation_config.{hpp,cpp} | 観測制約 + 重要度サンプリング範囲 |
| mass_function.{hpp,cpp} | レンズ質量サンプリング・残骸進化・等級変換 |
| event_sampler.{hpp,cpp} | MCループ・統計・stdout出力 |
| velocity_distribution.{hpp,cpp} | 速度サンプリング薄ラッパー |
| initialize.{hpp,cpp} | モデルオプション読み込み・RunContext初期化 |
| run_context.hpp | シミュレーション可変状態の集約 |
| internal/runtime.hpp | active_state宣言・互換マクロ (過渡的) |

### 互換ブリッジ (active_state)

`extern RunContext *active_state` + マクロ群が全 runtime .cpp ファイルで生きている。
各メソッドの先頭で `active_state = &ctx` をセットすることで科学計算部分の式は変えずに動く。
これは次の廃止ターゲット。

## 次のステップ

1. active_state マクロの廃止 — runtime 関数を明示的参照を取るクラスに変換
2. Python 直接パス — EventSampler が SimulationResult を直接返す
3. pre_gapmoe rewrite — gapmoe_src/ の書き直し (.gaplan/gapmoe_rewrite.md 参照)
