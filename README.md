---
name: dockq-benchmark
description: 计算严格的抗体-抗原分组 DockQ，用于评估抗体-抗原结构预测；链映射来自 CSV，支持 233 抗体 benchmark 和自定义数据集
---

# Skill: dockq-benchmark

## 功能

为预测的抗体/TCR–抗原复合物计算**严格的抗体-抗原界面 DockQ** 分数。

与朴素 `DockQ --short` 的关键区别：
- **只**评估 binder-靶点界面（抗体/TCR vs 抗原）
- **排除**抗体-抗体（重链-轻链）等内部界面——它们会虚增平均 DockQ
- 对多拷贝晶体结构、model 多拷贝输出和歧义链映射具有鲁棒性

---

## 两条可用流水线

| 流水线 | 适用场景 | 鲁棒性 |
|----------|------------|------------|
| **Robust**（`scripts/run_dockq_robust.py`） | 任何新数据集的默认选择，尤其是 HLA / TCR / 多拷贝晶体 | 处理多拷贝 native + model、β2m / 非界面抗原链、歧义映射 |
| 旧版分组（`run_dockq_ab_ag_grouped.py`） | 老的 SABDAB benchmark（单拷贝 native、干净的 Fab） | 简单的序列匹配映射；多拷贝时会失败 |

**默认用 robust。**

---

## 前置条件

| 项目 | 说明 |
|------|--------|
| Conda 环境 | `base`（DockQ 要求 numpy<2；scipy 用于距离矩阵） |
| DockQ | `dockq==2.1.3` |
| Biopython | 用 `Bio.PDB` 解析结构 |
| Robust runner | `scripts/run_dockq_robust.py`（本 skill 自带） |
| 辅助脚本 | `/data/corp/chenhao.zhang/ProtenixDesign/tools/data/PTX_performance_test/run_dockq*.py` |

```bash
source /data/corp/chenhao.zhang/miniforge3/etc/profile.d/conda.sh
conda activate base
```

---

## Robust DockQ 流水线

### 为什么需要 robust 流水线？

评估真实 PDB native（尤其是 HLA、TCR-pMHC、抗体-MHC）时，朴素分组 DockQ 会失败：

| 失效模式 | 原因 | Robust 的修复 |
|--------------|-------|-----------|
| `missing_mapping_X` | native 链数多于 model（多拷贝晶体） | 每个序列组选一条链 → 最小复合物 |
| `dockq_calc_failed`（无界面） | 抗原链在 native 中不接触抗体（如 β2m 在 HLA 另一侧） | 将抗原链裁剪为与抗体有足够界面的链（≥20 个原子接触） |
| 混合拷贝映射 | 纯按序列匹配时，不同 native 链会映射到不同 model 拷贝 | 以第一条抗体链为锚点 → 在 model 接触图中 BFS 圈定单个拷贝 |
| 抗原选择选错拷贝 | 晶体堆积使拷贝 2 的链"接触"拷贝 1 的锚点 | 用接触数（界面面积），而不是最小距离 |
| 截断链匹配不上 | 晶体拷贝的 N 端切割常不同 | 用基于最长公共子串的 `sequence_similarity` |

### 算法概览

```
输入: native.pdb, model.pdb, mapping_csv (id, ab_chains, native_chains)

1. 解析 CSV → 得到 ab_csv（binder 链 ID）、ag_csv（抗原链 ID）
   注意：多拷贝 native 的 CSV 可能列出所有拷贝的全部链

2. 按序列对 native 链分组（基于最长公共子串的聚类，阈值 0.65）
   例如 4 拷贝的 TCR α/β/HLA/β2m/肽 → 5 个组，每组 4 条链

3. 选择最小复合物（每组一条链）：
   - 锚点 = CSV 中第一条抗体链（权威拷贝）
   - 贪心生长：对剩下的每个组，挑选与已选链接触数最多
     的链（按界面面积打分）
   - 先处理抗体组（保证 binder 一致性），再处理抗原组

4. 裁剪抗原界面：丢弃与抗体侧原子接触 <20 的抗原链
   （滤掉 β2m、支架蛋白、远离界面的链）

5. 写出 NATIVE 子集：只保留选中的链，保证 DockQ 计算干净

6. 映射 NATIVE → MODEL 链：
   a) 计算序列相似度矩阵（基于最长公共子串）
   b) 若 model 链数 ≥ native 链数的 2 倍：把映射限制在单个 model 拷贝
      - 锚点：第一条 native 抗体链 → 最佳 model 匹配
      - 在 model 接触图中 BFS（自动调节 min_contacts 阈值，
        从 20 → 1500，直到 BFS 得到的链数 ≤ 目标大小）
   c) 按相似度降序贪心指派（差距 ≥0.05 视为无歧义）
   d) 空间步骤：歧义链按与已锚定映射的邻近度指派
   e) 强制补全剩余链

7. 运行分组 DockQ：
   - 构建 native_grouped.pdb（抗体→A，抗原→B）
   - 构建 model_grouped.pdb（映射后的链：抗体 model→A，抗原 model→B）
   - 可选：运行 DockQ 的 fix_numbering.pl 对齐残基编号
   - python -m DockQ model_grouped.pdb native_grouped.pdb --mapping AB:AB

8. 输出中的诊断信息：
   - native 最小界面距离（若 >5Å，说明 native 中无聚合物接触 = 无法评分）
   - 链映射及方法（seq/spatial/forced）和相似度
   - 实际选用的抗体/抗原链 ID（对比 CSV 列出的全并集）
```

### 运行方式

**批量模式**（CSV 驱动，多 case）：

```bash
python3 scripts/run_dockq_robust.py \
    --pred_dir /path/to/predictions \
    --native_dir /path/to/native_pdbs \
    --mapping_csv /path/to/mapping.csv \
    --output /path/to/output.csv \
    --n_cpu 32
```

**内联模式**（单 case，无需 CSV —— 适用于设计抗体或临时评估）：

```bash
python3 scripts/run_dockq_robust.py \
    --single_pred /path/to/model.pdb \
    --single_native /path/to/native.pdb \
    --ab_chains H,L \
    --ag_chains A \
    --case_id my_case \
    --output /path/to/output.csv
```

两种模式都会把每个 case 的 `chain_remap.json` 写到 `<output>.chain_remap/`
（或用 `--remap_dir <DIR>` 指定），供下游诊断使用。

输出 CSV 各列含义：

| 列 | 含义 |
|--------|---------|
| `case_id` | 样本标识符 |
| `ab_native_csv`, `ag_native_csv` | CSV 原始链列表（所有拷贝） |
| `ab_native_used`, `ag_native_used` | 实际用于 DockQ 的链（每组一条，经接触裁剪） |
| `ab_model`, `ag_model` | 映射到的 model 链 ID |
| `chain_mapping` | 每条链：native→model（相似度, 方法） |
| `n_native_chains`, `n_model_chains`, `n_assembly_units` | 诊断信息 |
| `min_iface_dist_native`, `min_iface_dist_model` | 界面有效性检查 |
| `dockq` | DockQ 分数（失败为 -1） |
| `note` | 失败原因或特殊处理说明 |

### 验证结果

在 HLA_family_v5 benchmark（487 个 case，抗体/TCR-HLA 晶体复合物）上：

| 指标 | 旧版分组 | Robust |
|--------|---------------|--------|
| 总算分样本 | 458 | 453 |
| 有效 DockQ | 381 (83.2%) | 445 (98.2%) |
| 平均 DockQ | 0.596 | **0.656** |
| 中位 DockQ | 0.717 | **0.807** |
| 高质量占比（≥0.80） | 38.6% | **51.5%** |
| 失败占比（<0.23） | 19.7% | 18.9% |

72 个 case 被修复（平均 DockQ 0.77），54 个提升 >0.1，23 个下降 >0.1，其余稳定。

---

## Mapping CSV 格式

```csv
id,ab_chains,native_chains
SABDAB00004,B,A
SABDAB00097,"C,D",E
SABDAB00492,"H,L",B
3PQY,"D,E,I,J,N,O,S,T","A,B,C,F,G,H,K,L,M,P,Q,R"   # 多拷贝：列出所有链
```

- `ab_chains`：native PDB 中 binder 侧链 ID（Fab 的 H+L、scFv、纳米抗体、TCR α+β）
- `native_chains`：靶点侧链 ID（HLA 场景：重链 + β2m + 肽；普通抗原：抗原链）
- 多拷贝晶体：列出**所有**拷贝的全部链；robust 流水线会自动选择一个最小复合物

### 自动构建 mapping CSV

#### 方案 A —— 从 RCSB 自动构建（已入库 PDB 的推荐做法）

对于任何已入库 RCSB 的 PDB ID，每个实体的描述（`pdbx_description`，
例如 `"2B04 heavy chain"`、`"Spike protein S1"`）和链 ID
（`entity_poly.pdbx_strand_id`）都可以直接从 Data API 获取：

- Entry 端点：`https://data.rcsb.org/rest/v1/core/entry/{PDB}` → `polymer_entity_ids` 列表
- 单实体端点：`https://data.rcsb.org/rest/v1/core/polymer_entity/{PDB}/{EID}`

对 `pdbx_description` 做关键词分类（heavy/light/Fab/VHH/nanobody/scFv/
immunoglobulin → 抗体；其余 → 抗原），再加上 SAbDab 注释检查，
远比序列前缀启发式可靠。已封装成现成脚本：

```bash
# 基本用法
python3 scripts/build_mapping_from_rcsb.py \
    7k9i 7lop 7str 6cw9 \
    --output mapping.csv --verbose

# 按 PDB 指定角色覆盖（设计的/非标准 binder，如 DARPin）
python3 scripts/build_mapping_from_rcsb.py 5o45 \
    --override "5O45:ab=B;ag=A" \
    --output mapping.csv

# 多个覆盖 + ID 文件
python3 scripts/build_mapping_from_rcsb.py --id-file pdb_ids.txt \
    --override "5O45:ab=B;ag=A" \
    --override "8U7M:ab=A;ag=B,C" \
    --output mapping.csv
```

**分类信号（按优先级）：**

1. **SAbDab 注释** —— 实体上的 `SABDAB_*` 特征（通过
   `rcsb_polymer_entity_feature_summary` 查询）。强阳性：该实体是抗体。
   比 RCSB 描述更新更快。
2. **描述关键词** —— heavy/light/Fab/VHH/scFv/nanobody/immunoglobulin
3. **TCR 描述** —— 无抗体时的兜底 binder（超出严格抗体-抗原流水线的
   范围；该行会被标记）。
4. 其余 → 抗原。

`--override` 是硬覆盖，绕过所有分类。用于不匹配标准抗体关键词
且尚未收录进 SAbDab 的设计 binder（DARPin、monobody、scFv-Fc 融合、
BiTE 等）。

输出列：`id, ab_chains, native_chains, ab_type, binder_kind, notes`。
前三列可直接喂给 DockQ runner；`ab_type`（VHH/Fab/TCR/Override）和
`binder_kind`（Ab/TCR/Override）是供过滤用的额外元数据。

**已验证 case（2026 年 4 月）：**

| PDB | 真实情况 | 脚本输出 | 状态 |
|-----|-------|---------------|--------|
| 7K9I | Fab(H,L) + Spike(A) | `ab=H,L ag=A type=Fab` | 正常 |
| 7STR | Fab(H,L) + NP(C) | `ab=H,L ag=C type=Fab` | 正常 |
| 7L5B | Fab(H,L) + Spike(A) | `ab=H,L ag=A type=Fab` | 正常 |
| 7OLZ | 2 个纳米抗体 + Spike | `notes=multiple_antibodies(H=2,L=0);needs_manual_split` | 已标记 |
| 7E9Q | 3 个 Fab + Spike 三聚体 | `notes=multiple_antibodies(H=3,L=3);needs_manual_split` | 已标记 |
| 7LOP | 2 个不同 Fab + 2 个 Spike | `notes=multiple_antibodies(H=4,L=4);needs_manual_split` | 已标记 |
| 8DF5 | 2 个 Fab + Spike + ACE2 | `notes=multiple_antibodies(H=4,L=4);needs_manual_split` | 已标记 |
| 6CW9 | TCR-CD1d（无抗体） | `ab=C,D ag=A,B type=TCR`（通过 `binder_kind=TCR` 标记，超出抗体-抗原范围） | 已标记 |
| 6VW1 | ACE2-RBD（无抗体） | `notes=no_binder` | 已标记 |
| 5O45 | PD-L1 + 大环肽（无抗体） | 需要 `--override "5O45:ab=B;ag=A"` | 手动 |

`notes` 非空的行在运行 DockQ 前**必须**人工检查：

- `no_antibody_chain[*]` —— 不是抗体-抗原复合物；移除该行，或改用
  其他评分方案。TCR-pMHC / 受体-配体等超出本 skill 范围
  （计划作为未来的子 skill）。
- `multiple_antibodies(...);needs_manual_split` —— 不对称单元中有
  2 个以上 Fab/VHH 拷贝。拆成每行一个 Fab，并为每行提供一份提取出的
  native PDB（用 `pdb-tools pdb_selchain`）。以 7LOP 为例：
  ```csv
  7LOP_CR3022,"W,V",Z
  7LOP_CV05,"X,Y",A
  ```
- `unknown_ab_subtype` —— 描述匹配了 "antibody/immunoglobulin" 但
  无法推断 H 还是 L；需人工核实。

#### 方案 B —— 设计抗体 / 非 RCSB 结构

对于没有入库 PDB 记录的设计，由用户在运行时提供链信息：

- **单 case**：直接把 `--ab_chains H,L --ag_chains A` 传给 runner
  （P1d）。完全绕过 CSV。
- **多 case**：按同样的 `id,ab_chains,native_chains` 格式写 CSV；
  设计方通常知道这些标识符。

流水线**不会**对设计序列做序列前缀启发式（`QVQL...` 等），因为设计
抗体经常不符合标准的胚系起始。HLA 家族 / TCR-pMHC 复合物超出本
skill 范围，相关旧版辅助脚本见
`/data/corp/chenhao.zhang/ProtenixDesign/input/HLA_family_v5/build_mapping.py`。

---

## 旧版：233 抗体 Benchmark 快速运行

对于标准的 233 抗体 Fab benchmark（单拷贝 native），用预置脚本：

```bash
python3 /data/corp/chenhao.zhang/protenix-ailux-v2/slurm/test_r2/run_dockq_all.py \
    --source sample1_seed30 --n_cpu 16
```

可用的 source：`sample1_seed1`、`sample1_seed5`、`sample1_seed10`、`sample1_seed15`、`sample1_seed30` 或 `all`。

| 项目 | 路径 |
|------|------|
| Mapping CSV | `/data/corp/chenhao.zhang/ProtenixDesign/tools/data/PTX_performance_test/mapping_233ab_full.csv` |
| Native PDB | `/data/corp/chenhao.zhang/ProtenixDesign/tools/data/PTX_performance_test/233antibody_native/` |

---

## 结果解读

### DockQ 质量分档

| 档位 | DockQ 范围 | 含义 |
|------|-------------|---------|
| **High** | ≥ 0.80 | 接近天然的对接 |
| **Medium** | 0.49 – 0.80 | 可接受的对接 |
| **Acceptable** | 0.23 – 0.49 | 勉强可用 |
| **Incorrect** | < 0.23 | 对接失败 |

### 失效模式（`note` 列）

| Note | 含义 | 可能的修复 |
|------|---------|--------------|
| `csv_chains_not_found_in_native` | CSV 中的链在 native PDB 里不存在 | 核实 CSV 链 ID 与 native 一致；链 ID 区分大小写 |
| `select_minimal_failed` | 无法选出抗体/抗原组（可能是 CSV 错了或没检测到分组） | 检查 ab_csv 至少有一条链存在于 native，ag 同理 |
| `missing_mapping_X` | native 链 X 无法映射到 model | model 缺了该实体；若是多拷贝问题，检查 `ab_native_used` |
| `native_no_contact(N.NA)` | native 中抗体和抗原链不接触（>5Å） | 可能是代谢物/辅因子桥接（如 MR1-RL 配体-MAIT TCR）—— DockQ 不适用 |
| `dockq_calc_failed` | DockQ 二进制失败；检查链编号、残基范围 | 尝试 DockQ 根目录下的 `fix_numbering.pl` |
| `relaxed_iface=N` | 界面裁剪用了 N-Å 阈值（>5Å） | 辅因子介导的相互作用 |
| `relaxed_contact=N` | 组装体提取用了 N-Å 阈值 | 松散复合物 |

### DockQ 根本不适用的情况

- **MR1-代谢物-MAIT TCR**：TCR 不接触 MR1 聚合物（>5Å），只通过小分子连接
- **纯 HLA-肽复合物**：没有抗体/TCR（按标题过滤的误报）
- **严重无序/切割的结构**：native 与 model 残基错配 >50%

这类情况 robust 流水线会报 `note=native_no_contact`，过滤掉即可。

---

## 内部机制小结

```
Native（多拷贝）                Model（多拷贝？）
┌──────────────────┐           ┌──────────────────┐
│ A B C D E (cp1)  │           │ A B C D E (cp1)  │
│ F G H I J (cp2)  │           │ F G H I J (cp2)  │
│ K L M N O (cp3)  │           └──────────────────┘
│ P Q R S T (cp4)  │
└──────────────────┘
        │
        ▼ select_minimal_complex（锚点 = CSV 中第一条抗体链）
        ▼ 将抗原裁剪到有足够的接触
        ▼
┌──────────────────┐
│ A B C D E (cp1)  │  ← 用于 DockQ 的 native 子集
└──────────────────┘
        │
        ▼ build_chain_mapping_with_spatial
        ▼   - 通过自动调参的 BFS 把 model 限制在单个拷贝
        ▼   - 序列 + 空间指派
        ▼
   {A→A, B→B, C→C, D→D, E→E}（在 model 单拷贝内一致）
        │
        ▼ build_grouped_pdb（抗体→A，抗原→B）
        ▼ DockQ --mapping AB:AB
        ▼
     DockQ 分数
```

---

## 局限性 & 结果不可信的情况

本流水线给出的抗体-抗原界面 DockQ 是诚实的，但有几个结构性边界
情况，解读数字前请注意：

### 1. 同源多聚体抗原 + 单个 binder（如 spike 三聚体 + 1 个 Fab）

当抗原是同源多聚体（多条相同序列的链，如 spike 三聚体、GPCR 二聚体、
病毒衣壳）时，只有**一个**拷贝真正接触 binder。选哪个拷贝并不平凡：

- **Native 侧**：流水线挑选与锚定抗体接触原子数最多的抗原链
  —— 这是正确的。
- **Model 侧（2026-04 的 P1 修复）**：抗原链在 model 中按同样的
  "与 model 抗体接触最多"规则独立重选，而不是按与 native 锚定
  抗体映射的空间邻近度。这样：
  - 对接正确 → model 选中同一条逻辑抗原链 → DockQ 诚实。
  - 对接错误（binder 对到非天然界面或错误拷贝）→ model 选中它
    实际接触的那条抗原链；若一条都不接触，报 `note=model_no_iface`
    且 DockQ ≈ 0。这是**真实**的失败信号。
- **不要**枚举所有抗原链取最大 DockQ —— 当 model 对接错误时，
  那种策略可能伪映射到一条未被接触的抗原链，报出空界面
  （DockQ 未定义 / 被残端原子虚增）。

### 2. 一个不对称单元中有多拷贝抗体（如 7LOP、8DF5、7E9Q）

如果 native PDB 含有 2 个以上结合在抗原上的 Fab/VHH 拷贝——无论
它们是**同一抗体**的多拷贝（如 7E9Q：spike 三聚体上 3 个相同 Fab）
还是**不同抗体**（如 7LOP、8DF5）——mapping CSV 必须拆成每行一个
Fab。每行需要自己的 native PDB（用 `pdb-tools pdb_selchain` 提取）。
自动构建脚本（`build_mapping_from_rcsb.py`）会用
`notes=multiple_antibodies(H=N,L=N);needs_manual_split` 标记这类情况
—— 不要把带该标记的行直接喂给 runner。

**为什么同一抗体多拷贝也必须拆分**：runner 的 `select_minimal_complex`
以 CSV 中第一条抗体链为锚点，再按接触数贪心加链；在同源多 Fab 情况
下，这可能把拷贝 1 的 L 链和拷贝 2 的 H 链配成一对（在 7E9Q 的
self-vs-self 上已验证会失败 → DockQ=0.32）。预先拆分可消除该歧义。

### 3. RCSB 关键词分类器有长尾盲区

`build_mapping_from_rcsb.py` 对 `pdbx_description` 做关键词匹配
（heavy/light/Fab/VHH/nanobody/scFv/immunoglobulin）。它能正确标记
7K9I / 7STR / 7LOP / 8DF5 / 6CW9 / 6VW1（已验证）。但以下描述不匹配
标准抗体关键词，**会被误分类为抗原**：

- `BiTE`（双特异性 T 细胞衔接器）
- `DARPin`、`affibody`、`monobody`、`anticalin`
- `engineered T-cell receptor`
- 某些嵌合/融合构建体（如 "scFv-Fc fusion to receptor X"）

**缓解措施**：脚本还会查询 `rcsb_polymer_entity_annotation` 中的
`SAbDab` 类型注释（P1b 修复）；SAbDab 覆盖了大多数工程化抗体形式。
仍不确定的一律用 `--override` 参数（P1c）或手动编辑 CSV。

### 4. 设计抗体 / 非 RCSB 结构

自动映射流水线只支持 RCSB。对于没有入库 PDB 记录的设计抗体：
- 运行时通过 `--ab_chains H,L --ag_chains A` 提供链 ID（P1d，单 case）
- 或提供符合 `id,ab_chains,native_chains` 格式的 CSV。

流水线**不会**对预测结构做序列前缀启发式（`QVQL...`），因为设计
抗体经常不符合标准的胚系起始。

### 5. 分组链 ID 在 DockQ 二进制输出中丢失原始映射

流水线内部会把抗体重命名为 A、抗原重命名为 B，用于
`DockQ --mapping AB:AB`。因此 DockQ stdout 中的逐残基/逐界面残基
诊断信息引用的是 A/B，而不是原始链 ID。完整的原始链映射写在
`<output>.chain_remap/<case_id>.json`（P2 修复）供诊断。

---

## 常用操作建议

1. **务必核实 CSV 链 ID 与 native PDB 一致**（区分大小写；多字符链 ID 必要时通过 `cif→pdb` 转成单字母）
2. **超大数据集**（>1000 case）用 `--n_cpu` 匹配集群核数；H100 节点上约 10 秒/样本
3. **检查低 DockQ case 的 `chain_mapping` 和 `*_used` 列** —— 链选错是假性低分的头号来源
4. **β2m 的位置很重要**：如果 CSV 把 β2m 放在抗体侧，robust 流水线可能锚错。一律把 β2m 归到抗原侧
5. **未知数据集**：先用默认 `contact_thresh=5.0` 和 `min_contact_atoms=20`；只有当很多 case 出现 `relaxed_iface` 标记时才调参
