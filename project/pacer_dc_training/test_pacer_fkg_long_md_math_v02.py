import numpy as np

from pacer_fkg_long_md_math_v02 import (
    AXIS_SIGNS,
    preprocess_features,
    rff_map,
    temporal_block_mean,
    signed_block_contrast,
    region_block_vectors,
    pooled_replica_vector,
    summarize_three_replicas,
)


def assert_close(a, b, atol=1e-12):
    np.testing.assert_allclose(a, b, atol=atol, rtol=0)


# ------------------------------------------------------------
# 1. Frozen sign definitions
# ------------------------------------------------------------
assert AXIS_SIGNS["dAGO"] == {
    "candidate": 1.0,
    "apo": -1.0,
}
assert AXIS_SIGNS["dPAM"] == {
    "candidate_probe": 1.0,
    "probe": -1.0,
}
assert AXIS_SIGNS["dINT"] == {
    "candidate_probe": 1.0,
    "probe": -1.0,
    "candidate": -1.0,
    "apo": 1.0,
}


# ------------------------------------------------------------
# 2. Exact synthetic signed contrasts
#
# apo = 1
# candidate = 4      => dAGO = +3
#
# probe = 10
# candidate_probe=17 => dPAM = +7
#
# dINT = dPAM - dAGO = +4
# ------------------------------------------------------------
shape = (5, 3, 2)

ctx = {
    "apo": np.full(shape, 1.0, dtype=np.float32),
    "probe": np.full(shape, 10.0, dtype=np.float32),
    "candidate": np.full(shape, 4.0, dtype=np.float32),
    "candidate_probe": np.full(shape, 17.0, dtype=np.float32),
}

ago = signed_block_contrast(ctx, "dAGO")
pam = signed_block_contrast(ctx, "dPAM")
inter = signed_block_contrast(ctx, "dINT")

assert_close(ago, 3.0)
assert_close(pam, 7.0)
assert_close(inter, 4.0)
assert_close(inter, pam - ago)


# ------------------------------------------------------------
# 3. Paired-block cancellation
#
# Add the same arbitrary block-dependent common background
# to all four contexts. All three signed contrasts must remain
# unchanged because each contrast is evaluated at the same
# paired block index.
# ------------------------------------------------------------
rng = np.random.default_rng(12345)
background = rng.normal(size=shape).astype(np.float32)

ctx_bg = {
    k: v + background
    for k, v in ctx.items()
}

assert_close(
    signed_block_contrast(ctx_bg, "dAGO"),
    ago,
    atol=2e-6,
)
assert_close(
    signed_block_contrast(ctx_bg, "dPAM"),
    pam,
    atol=2e-6,
)
assert_close(
    signed_block_contrast(ctx_bg, "dINT"),
    inter,
    atol=4e-6,
)


# ------------------------------------------------------------
# 4. Pairing is index-sensitive.
# Permuting one context must alter a time-varying paired
# contrast; the code must NOT silently sort/re-pair blocks.
# ------------------------------------------------------------
block_signal = np.arange(5, dtype=np.float32)[:, None, None]
dynamic = {
    "apo": np.zeros(shape, dtype=np.float32),
    "probe": np.zeros(shape, dtype=np.float32),
    "candidate": np.broadcast_to(
        block_signal,
        shape,
    ).copy(),
    "candidate_probe": np.broadcast_to(
        2 * block_signal,
        shape,
    ).copy(),
}

correct = signed_block_contrast(dynamic, "dINT")

mispaired = dict(dynamic)
mispaired["candidate_probe"] = np.roll(
    dynamic["candidate_probe"],
    1,
    axis=0,
)

wrong = signed_block_contrast(mispaired, "dINT")

assert not np.array_equal(correct, wrong)


# ------------------------------------------------------------
# 5. Region averaging and temporal pooling
# ------------------------------------------------------------
nodes = np.zeros((5, 4, 2), dtype=np.float32)
nodes[:, 1, :] = [2.0, 4.0]
nodes[:, 3, :] = [6.0, 8.0]

regional = region_block_vectors(nodes, [1, 3])
assert_close(
    regional,
    np.tile([4.0, 6.0], (5, 1)),
)

pooled = pooled_replica_vector(regional)
assert_close(pooled, [4.0, 6.0])


# ------------------------------------------------------------
# 6. Three-replica direction summary.
# No qualification gate is inferred.
# ------------------------------------------------------------
summary = summarize_three_replicas({
    1: np.array([1.0, 0.0]),
    2: np.array([2.0, 0.0]),
    3: np.array([3.0, 0.0]),
})

assert_close(
    summary["pairwise_direction_cosines"]["R1_vs_R2"],
    1.0,
)
assert_close(
    summary["pairwise_direction_cosines"]["R1_vs_R3"],
    1.0,
)
assert_close(
    summary["pairwise_direction_cosines"]["R2_vs_R3"],
    1.0,
)
assert summary["qualification_gate"] == "NOT_DEFINED"

opposed = summarize_three_replicas({
    1: np.array([1.0, 0.0]),
    2: np.array([2.0, 0.0]),
    3: np.array([-1.0, 0.0]),
})

assert_close(
    opposed["pairwise_direction_cosines"]["R1_vs_R3"],
    -1.0,
)
assert opposed["qualification_gate"] == "NOT_DEFINED"


# ------------------------------------------------------------
# 7. Frozen preprocessing semantics:
# per-frame residue centering happens BEFORE med/scale.
# ------------------------------------------------------------
raw = np.array(
    [
        [[1.0, 10.0],
         [3.0, 14.0]],
    ],
    dtype=np.float32,
)
med = np.array([0.5, -1.0], dtype=np.float64)
scale = np.array([2.0, 4.0], dtype=np.float64)

pre = preprocess_features(raw, med, scale)

centered = raw - raw.mean(axis=1, keepdims=True)
expected = ((centered - med) / scale).astype(np.float32)

assert_close(pre, expected)


# ------------------------------------------------------------
# 8. Generic RFF formula
# ------------------------------------------------------------
x = np.array([[[1.0, 2.0]]], dtype=np.float32)
w = np.array(
    [[1.0, 0.0],
     [0.0, 1.0]],
    dtype=np.float32,
)
b = np.array([0.0, 0.0], dtype=np.float32)

z = rff_map(x, w, b)

expected_z = np.sqrt(2.0 / 2.0) * np.cos(
    np.array([1.0, 2.0], dtype=np.float32)
)

assert_close(z[0, 0], expected_z, atol=1e-7)


# ------------------------------------------------------------
# 9. Temporal block partition preserves order
# ------------------------------------------------------------
x = np.arange(
    40,
    dtype=np.float32,
).reshape(10, 2, 2)

blocked = temporal_block_mean(x, block_frames=2)

assert blocked.shape == (5, 2, 2)
assert_close(blocked[0], x[0:2].mean(axis=0))
assert_close(blocked[4], x[8:10].mean(axis=0))


print("PACER_FKG_LONG_MD_SYNTHETIC_PASS")
print("signs_dAGO", AXIS_SIGNS["dAGO"])
print("signs_dPAM", AXIS_SIGNS["dPAM"])
print("signs_dINT", AXIS_SIGNS["dINT"])
print("paired_block_indexing PASS")
print("common_background_cancellation PASS")
print("region_pooling PASS")
print("three_replica_aggregation PASS")
print("opposed_replica_cosine PASS")
print("preprocessing_order PASS")
print("rff_formula PASS")
print("qualification_gate NOT_DEFINED")
print("kernel_u2_used False")
