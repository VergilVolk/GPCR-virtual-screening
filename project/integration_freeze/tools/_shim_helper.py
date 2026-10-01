from pathlib import Path
import sys
REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))


class _Shim:
    n_frames = 400
    block_frames = 20
    blocks_per_trajectory = 20
