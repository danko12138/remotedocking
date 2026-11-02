# ABCB1 / ABCG2 — ZINC22 remote docking

利用 Churchill 现有 ZINC22 分片筛选 ABCB1、ABCG2 候选。实验活性标签用于方法校验和解释，**不作为分子入队条件**。本仓库是可运行的小批验证起点，提供准备、Vina CPU / GNINA CPU-GPU 入口、逐任务断点记录和结果表；尚未对 Churchill 验证吞吐量，不能把它当作已验收的百万级并行系统。

## 已带入的材料

- ABCB1：9CR8、9CTC、9CTF；ABCG2：6VXF、6VXJ、7OJ8。
- 每套受体的原有 PDBQT、搜索盒、SHA256、坐标框架和历史限制。
- 十个现成参考输入：tariquidar、zosuquidar、elacridar、KO143、mitoxantrone、topotecan、erlotinib、gefitinib、imatinib、sulfasalazine。具体盐型/立体身份以 `data/references/manifest.json` 完整键为准，名称不用于合并身份。
- 参考输入的已有准备收据。部分源结构立体化学未完全指定，代表构象不等于实验赋值。

现有受体包含 apo 结构、原生姿势恢复未通过项，以及未完整验证的质子化约定，详情见 `data/receptors/manifest.json`。它们支持工程衔接，**不代表新筛选方法已获得科学验证**。当前没有承诺原生配体重对接 RMSD；原始 PDB/CIF、共晶配体与对应残基映射仍需补充到方法验证包。

## Churchill 上的起步步骤

```bash
git clone git@github.com:danko12138/remotedocking.git
cd remotedocking
conda env create -f environment.yml
conda activate abc-zinc22
python -m pip install --no-build-isolation --no-deps -e .
cp configs/churchill.example.json configs/local.json
```

没有 SSH key 时可用 HTTPS 克隆：`git clone https://github.com/danko12138/remotedocking.git`。

**没有 Conda，或下载报 `Could not resolve host: conda.anaconda.org`：** 可用独立安装脚本，从官方 GitHub Release 下载固定版本 Micromamba 并验证 SHA256，在个人目录安装环境。选择清华 conda-forge 镜像绕过该 Anaconda 域名：

```bash
git pull --ff-only
python3 tools/setup_environment.py --channel tuna
```

默认安装到仓库上级目录的 `.abc-zinc-tools` 与 `envs/abc-zinc22`；不修改 DNS、`.bashrc`、系统 Python 或共享程序，也不启动 docking。清华镜像的可达性和版本同步情况仍需远程实际确认；脚本遇到失败即停止，保留错误。可先加 `--dry-run` 查看执行计划。包版本不会因换镜像而自动放宽。

安装后无需 activate，可使用：

```bash
../.abc-zinc-tools/bin/micromamba-2.10.0-0 --no-rc \
  -r ../.abc-zinc-tools/root run -p ../envs/abc-zinc22 zincscreen --help
```

来源：[Micromamba 官方 GitHub Releases](https://github.com/mamba-org/micromamba-releases/releases/tag/2.10.0-0)、[清华镜像配置说明](https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/)。

`environment.yml` 固定了本地测试所用的主要版本，**尚未在 Churchill 求解安装，不是远程平台锁文件**。可使用现有兼容环境，先执行自检并留存实际包版本。GNINA 为外部程序，仓库和环境文件没有携带它的二进制或 CUDA 库。

编辑 `configs/local.json`：

| 字段 | Churchill 设置 |
| --- | --- |
| `workspace` | `/home/zyc/workspaces/member26/abcb1_abcg2_zinc22_work`，独立保存任务与断点 |
| `zinc_root` | `/home/zyc/projects/uhrf1_usp7/data/compounds/zinc22_broad/raw`，只读原始库 |
| `engine.executable` | 真实 GNINA 可执行文件的绝对路径，或 PATH 中的命令 |
| `engine.backend` / `use_gpu` | GNINA GPU：`gnina` / `true`；Vina CPU 基线：`vina` / `false` |
| `engine.cpu_threads` | 初始 2；结合共享机器实际分配调整 |
| `engine.device` | `nvidia-smi` 对应的物理 GPU 序号，初始 0 |
| `gpu` | 提交前空闲显存和利用率阈值，按共享使用安排调整 |

资源自检（不提交 docking）：

```bash
lscpu
free -h
nvidia-smi
df -h /home/zyc/workspaces/member26
zincscreen doctor --config configs/local.json --output work/doctor.json
```

自检失败返回退出码 2，逐项原因保存在报告中。示例默认选择 GNINA GPU；如果暂时只有环境内的 Vina，改为 `backend=vina`、`executable=vina`、`use_gpu=false` 后先验证 CPU 基线。不能给 Vina CPU 配置 `use_gpu=true`。

用户硬件快照：WSL 可见 20 CPU、94 GiB RAM、GPU 97,887 MiB、工作盘剩余约 1.1 TiB；GPU 当时 **89% 利用率且已有多项 GNINA**。默认配置会等待 GPU 低于阈值，并保留 50 GiB 磁盘余量。GNINA 完整路径和远程依赖仍需操作者确认，详见远程交接文档。

## 第一轮参考配体与小批测试

```bash
# 复用已带入的十个参考输入，重新检查文件与原子映射；不重新生成参考构象。
zincscreen prepare --config configs/local.json --references

# 默认串行、20 项：10 配体 × 2 受体；此处仅作工程验证。
zincscreen run --config configs/local.json \
  --prepared /home/zyc/workspaces/member26/abcb1_abcg2_zinc22_work/prepared/references.jsonl \
  --receptors 9CTF 6VXJ --seed 42 --max-jobs 20
```

扫描文件清单和有界采样：

```bash
zincscreen catalog --config configs/local.json --output work/zinc_catalog.jsonl
zincscreen sample --config configs/local.json --output work/zinc100.jsonl \
  --limit 100 --per-file 4 --scan-lines 10000 --seed 42
zincscreen prepare --config configs/local.json --input work/zinc100.jsonl
```

准备命令最后输出 `manifest` 路径。将该路径作为 `--prepared`，先运行 20 项，再扩大 `--max-jobs` 或选全部六受体：

```bash
zincscreen run --config configs/local.json --prepared <manifest路径> \
  --receptors 9CR8 9CTC 9CTF 6VXF 6VXJ 7OJ8 --seed 42 --max-jobs 600
zincscreen report --config configs/local.json --output work/scores.csv
```

`sample` 在可识别的 Hxx 路径分层间轮流访问，每文件读取至多指定行数，在该前缀中做 reservoir sampling。**这是有界工程样本，不是全库均匀样本或完整分子数审计**；收据保留访问文件、文件 stat、抽样参数与异常示例。`catalog` 只统计文件和压缩字节，历史 1,383 文件 / 1.29 TiB 尚待实际对账。文件后缀可读取 gzip/bzip2/xz/未压缩 SMILES；首轮需要核实库的实际格式。ZINC ID 与原始 SMILES 一起保留，不从活性标签补造身份。

## 恢复、更新和数据保留

- 相同输入、程序哈希、参数、seed 和代码版本的已完成任务，检查输出哈希后复用。
- 新 seed/新程序/新盒/代码改变形成不同任务身份，不覆盖旧结果；因此升级代码后不要误把新任务当成原任务续跑。
- 失败默认保留，不循环重试。加 `--retry-failed` 明确处理失败项；`ANALYSIS_FAILED` 仅重解析已留存输出。
- 被中断或遗留的 `RUNNING` 记录也默认保留。确认旧计算确实已停止后再用 `--retry-failed`；不得在不同机器共用同一个 workspace。
- 单进程文件锁保护本机同一 workspace。默认串行，不创建定时器、后台服务或自动开机任务。
- GPU 提交前检查显存及利用率，不终止其他用户进程。瞬时检查不是资源预留，需遵守工作站共享安排；使用本入口时清除外部 `CUDA_VISIBLE_DEVICES`，以配置物理序号为准。
- ZINC、运行结果、环境、缓存及 `configs/local*.json` 不进入 Git。`git pull` 只更新代码与输入包，工作目录独立保存。

## 输出与筛选口径

每个任务保存 `SPEC.json`、每次独立尝试、实际命令、程序/输入哈希、日志、原始姿势、全部已解析能量、几何标记及近邻残基。表格分别记录最低原始分和最低几何合格姿势分。正能量不被隐藏，也不强制改成负值；零个合格姿势不因极低分而放行。

GNINA 的经验能量、CNNscore、CNNaffinity 分开保存。4 Å 内残基列表只是距离接触，不能叫氢键或实验关键相互作用。当前 QC 是粗筛：有限坐标、原子/结构映射、键长和严重受体碰撞；没有完整完成非键合内部碰撞、所有立体构型、共价化学或氢键验证，仍需精细检查与可视化。

本轮仍是刚性受体＋可动配体搜索。多构象不是连续蛋白动力学；AF3/MD 未接入本入口。跨靶点/后端原始分数不可直接标为实验选择性，不将历史 UHRF1–USP7 A/B 分数用于本项目排名。

后续顺序：参考口袋验证 → 分层 1,000 分子六构象测速 → 按实测扩大筛库 → 更充分搜索/多 seed/姿势与化学检查 → 多样化候选及独立复筛。百万级分片导入、跨主机调度、Uni-Dock/Vina-GPU/AD4 适配和 AF3/MD 属于后续阶段，本版未实现。

## 测试与参考

```bash
python -m pytest -q
```

[远程执行步骤与硬件信息](docs/CHURCHILL.md) · [验证记录](docs/VALIDATION.json)

接口依据：[AutoDock Vina](https://autodock-vina.readthedocs.io/en/latest/docking_basic.html)、[GNINA](https://github.com/gnina/gnina)、[RDKit Gasteiger](https://www.rdkit.org/docs/source/rdkit.Chem.rdPartialCharges.html)。化学准备沿用本地项目的身份/立体检查及 OpenBabel 转换思路，新增输入要求力场支持且优化收敛；电荷使用 RDKit 严格 Gasteiger（缺参报错），核对原子后显式写入 OpenBabel，并检查 PDBQT 总电荷。第三方程序与库按其各自许可安装；仓库未打包程序、模型权重、实验标签全表或 ZINC 原库。
