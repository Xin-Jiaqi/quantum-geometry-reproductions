#!/bin/bash
#=============================================================================
# BD36 Phase 0 | para 就地接口验证单作业 (2253_bilayer_baseline)
#   R1.3 版：VASP2WAN90_v2_fix 补丁约束已核实 ——
#     * LWRITE_SPN 仅串行 (WDES%COMM%NCPU==1)，否则被静默跳过
#     * LSPN_FMTED=.TRUE. -> ASCII .spn；默认 binary (unformatted)
#     * 预处理(SOC 静态)与接口均 ISYM=-1 -> 同一完整 81 k 全 BZ 列表
#     * LWRITE_UNK=.FALSE. -> UNK 不是必需产物
#   资源上限: 28 核 x <=2 h = <=56 核时（单作业，Phase 0）
#   用法: 批准后, 在本目录 (INPUTS/2253_bilayer_baseline/) 执行
#         sbatch phase0_step0.sbatch
#   验收: 脚本末尾 ACCEPTANCE 汇总; 任一 EXIT-N 退出条件命中即停。
#=============================================================================
#SBATCH --job-name=bd36_p0_2253
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=root.bingxing
#SBATCH --output=phase0_%j.log

set -euo pipefail

STD="$(cd "$(dirname "$0")" && pwd)"          # structure dir (this INPUTS dir)
ROOT="$(cd "$STD/../.." && pwd)"              # research root
SCRIPTS="${SCRIPTS_DIR:-$ROOT/scripts}"
W90KPT_SCRIPT="${W90KPT_SCRIPT:-$SCRIPTS/make_win_kpoints.py}"
READSPN_SCRIPT="${READSPN_SCRIPT:-$SCRIPTS/read_spn.py}"
PYTHON="${PYTHON:-python3}"

# --- 环境装配 (run_env.sh 本地生成, 不入库; 提供 $VASP $W90 与 module load) ---
if [ -f "$STD/run_env.sh" ]; then source "$STD/run_env.sh"; fi
if [ -z "${VASP:-}" ] || [ -z "${W90:-}" ]; then
  echo "FATAL: run_env.sh must define VASP and W90 (see run_env.example.sh)"; exit 1
fi
export OMP_NUM_THREADS=1

# Fail closed on re-submission. Never re-use old .spn/.eig/.amn/.mmn
# from a previous attempt: each Slurm job gets an isolated directory.
RUN_DIR="$STD/phase0_run_${SLURM_JOB_ID:-manual}"
[ ! -e "$RUN_DIR" ] || { echo "FATAL: run directory already exists: $RUN_DIR (no duplicate run)"; exit 1; }
mkdir "$RUN_DIR" && cd "$RUN_DIR"
echo "== [P0.a] SOC 静态 (ISYM=-1, 9x9x1 -> 81 k 全 BZ; CHGCAR 供接口) =="
cp "$STD/INCAR_soc" INCAR
cp "$STD/KPOINTS" KPOINTS
cp "$STD/POSCAR" POSCAR
# Stage private POTCAR into the isolated run directory. Do not commit this file.
# run_env.sh may set POTCAR_SOURCE; alternatively place a private POTCAR beside this script.
if [ -n "${POTCAR_SOURCE:-}" ] && [ -f "$POTCAR_SOURCE" ]; then
  cp "$POTCAR_SOURCE" POTCAR
elif [ -f "$STD/POTCAR" ]; then
  cp "$STD/POTCAR" POTCAR
fi
[ -s POTCAR ] || { echo "FATAL: POTCAR missing; set POTCAR_SOURCE in private run_env.sh"; exit 1; }
mpirun -np 28 "$VASP" > soc_static.out 2>&1 || { echo "EXIT-1: SOC static failed"; tail -30 soc_static.out; exit 1; }
grep -q "k-points in reciprocal lattice" OUTCAR || { echo "EXIT-2: no k-points in OUTCAR"; exit 1; }

echo "== [P0.b] 提取实际 k 列表 -> win kpoints 块 (必须 81, Gamma 首点, 权重和=81) =="
[ -s IBZKPT ] || { echo "EXIT-3: IBZKPT not generated"; exit 1; }
"$PYTHON" "$W90KPT_SCRIPT" IBZKPT "9 9 1" > kpoints_block.txt || { echo "EXIT-3: kpoint extraction failed"; exit 1; }
grep -q "count=81" kpoints_block.txt || { echo "EXIT-4: k-count != 81"; exit 1; }

echo "== [P0.c] 组装 win + wannier90.x -pp =="
cp "$STD/wannier90.win" wannier90.win
"$PYTHON" - <<'EOF'
from pathlib import Path
import re
text = Path("wannier90.win").read_text()
text = re.sub(r"(?s)begin kpoints.*?end kpoints\n?", "", text)
block = Path("kpoints_block.txt").read_text()
Path("wannier90.win").write_text(text + "\n" + block)
n = sum(1 for l in block.splitlines() if len(l.split()) == 3)
print(f"win assembled with {n} k-lines")
EOF
"$W90" -pp wannier90 > w90_pp.out 2>&1 || { echo "EXIT-5: wannier90.x -pp failed (MPI/MKL env?)"; tail -20 w90_pp.out; exit 1; }

echo "== [P0.d] 接口作业 (LWANNIER90+LWRITE_SPN+LSPN_FMTED, ISYM=-1; 必须串行!) =="
cp "$STD/INCAR_soc_w90_vasp544" INCAR
cp "$STD/POSCAR" POSCAR; cp "$STD/KPOINTS" KPOINTS
[ -f CHGCAR ] || { echo "FATAL: CHGCAR missing after static"; exit 1; }
# WRITE_SPN 仅串行: 单 MPI rank 运行接口步
mpirun -np 1 "$VASP" > w90_iface.out 2>&1 || { echo "EXIT-6: interface run failed"; tail -40 w90_iface.out; exit 1; }
# Do not depend on a patch-specific stdout phrase; validate the output artifact.
[ -s wannier90.spn ] || { echo "EXIT-7: .spn missing (serial WRITE_SPN skipped or interface failed)"; grep -i spn w90_iface.out | tail -10 || true; exit 1; }

echo "== [P0.e] 真实验收 =="
# 必需产物: .eig .amn .mmn .spn (UNK 非必需, LWRITE_UNK=F)
for f in wannier90.eig wannier90.amn wannier90.mmn wannier90.spn; do
  [ -s "$f" ] || { echo "EXIT-8: required product $f missing/empty"; exit 1; }
done
# .spn 解析 (自动检测 ASCII/binary): 维数 (81,64,64,3), Hermite
"$PYTHON" "$READSPN_SCRIPT" wannier90.spn 81 > spn_check.json || { echo "EXIT-9: read_spn rejected .spn"; exit 1; }
"$PYTHON" - <<'EOF'
import json
r = json.load(open("spn_check.json"))
assert r["n_k"] == 81 and r["num_bands"] == 64, r
assert r["dims"] == [81, 64, 64, 3], r["dims"]
assert r["hermitian"] is True, "non-Hermitian .spn"
print("spn ACCEPTED:", r)
EOF
# k 一致性: 接口实际 IBZKPT 的完整 k 列表与 win kpoints 块一致
[ -s IBZKPT ] || { echo "EXIT-10: interface IBZKPT missing"; exit 1; }
"$PYTHON" "$W90KPT_SCRIPT" IBZKPT "9 9 1" > iface_k.txt
diff <(grep -v "^begin\|^end\|^#" kpoints_block.txt) <(grep -v "^begin\|^end\|^#" iface_k.txt) >/dev/null \
  || { echo "EXIT-10: interface k-list != win kpoints"; exit 1; }

echo "=========== PHASE 0 ACCEPTANCE: PASS ==========="
echo "products: wannier90.eig/.amn/.mmn/.spn; spn dims (81,64,64,3) Hermitian; k=81 consistent"
echo "next: Round-1 (Phase 1, 448 core-h) per PARA_STEP0_SCHEME.md / COST_AND_PERMISSION_GATE.md"
