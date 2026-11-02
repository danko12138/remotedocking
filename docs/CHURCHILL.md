# Churchill 交接

## 已知与待核实

用户提供：Ubuntu/WSL2 共享工作站、个人空间 `/home/zyc/workspaces/member26`；ZINC 根目录见示例配置。本轮没有成功登录或提交远程任务，以下为用户提供的终端快照，未由本代理直接复测。

| 项目 | 用户提供的硬件输出 |
| --- | --- |
| CPU | Core Ultra 9 285K，WSL 可见 20 个 CPU / 20 核；这不是宿主机物理核心总数审计 |
| 内存 | 94 GiB，总快照可用约 90 GiB；swap 32 GiB |
| GPU | RTX PRO 6000 Blackwell，97,887 MiB 显存，快照使用 3,325 MiB |
| GPU 活跃程度 | **89%**，已有 4 个 GNINA 与 1 个 Python 计算进程；不应因显存空闲而认定 GPU 空闲 |
| 驱动 | nvidia-smi 显示 driver 597.06，CUDA compatibility 13.2；不是已安装 CUDA toolkit 版本证明 |
| 工作盘 | 3.4 TiB，总已用 67%，可用约 1.1 TiB，挂载在 `/` |

主机快照时间字符串是 `Mon Nov 2 10:44:12 2026`，未提供时区，不将它转换成用户纽约时间或当作持续状态。进程列表出现 `gnina.cuda12.8.static`，说明该快照已有相关程序在运行，但没有给出可供本仓库使用的完整可执行文件路径、版本或环境。

按这份配置，首版继续 **1 个 docking 任务 × 2 个 CPU 线程**，准备单进程、4 GiB 限额。保留 50 GiB 磁盘余量；GPU 仅在利用率 ≤30% 且空闲显存 ≥12 GiB 时提交下一项。89% 快照下自检会报告 GPU busy，运行入口暂缓提交，不杀现有任务。这些是试跑初值，并非调度器给出的独占额度。

硬件信息由操作者执行：

```bash
lscpu
free -h
nvidia-smi
df -h /home/zyc/workspaces/member26
command -v conda micromamba gnina vina unidock
```

CPU 核心数用于决定 CPU/准备资源；可用 RAM 限制并发准备和 GNINA；GPU 显存与实际任务峰值决定 GPU 并发；磁盘空间和吞吐限制分片大小。96 GiB 不能直接换算同时运行数量。本版串行用于建立基线，收到硬件输出及实测后再扩展并发。

现在还需确认实际软件路径，可执行 `command -v gnina gnina.cuda12.8.static vina conda micromamba`。选择操作者有权使用的程序，把真实路径填入 `engine.executable`，再运行 doctor；不要仅根据进程名称猜测程序位于 `/`。

## 首批执行及反馈

1. clone，安装或选择已有环境，填写 `configs/local.json`。
2. 执行 `doctor`；保留 `conda list --explicit > work/environment-explicit.txt`。这是远程求解后形成的真实环境记录，不应把本地环境路径复制过去。
3. 参考配体 20 项小试；核对盒和接触、原始姿势，不把示例受体当作已验收的最佳口袋。
4. `catalog` 对账文件，再抽 100 个 ZINC 分子；检查准备失败原因，不以实验标签排除无标签分子。
5. 汇总端到端墙钟时间、引擎耗时、输入准备耗时、成功率、显存峰值及输出大小。外部 GPU 工具的采样窗口应注明；一次 `nvidia-smi` 不能代表全过程峰值。
6. 决定 1,000 分子六构象测试及后续扩展。先比较两个靶点各自结果，选择性作为后续经校准的比较。

## 边界

只在个人空间写输出；原始 ZINC 和其他项目库只读使用。不要将 UHRF1/USP7 的受体、搜索盒、A/B 分数沿用到 ABCB1/ABCG2。脚本不会自动创建云端任务、监控或持续下载数据。

仓库上传不需要 Churchill 的登录信息。GitHub SSH 克隆是否可用由 Churchill 的 GitHub key 决定，HTTPS 公有仓库可直接克隆；若为私有仓库，需要操作者自己的授权。没有要求将本机 SSH 私钥复制过去。
