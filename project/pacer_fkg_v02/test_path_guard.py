import tempfile
import unittest
from pathlib import Path

from path_guard import ALLOWED_RELATIVE_ROOTS, PROTECTED_RELATIVE_ROOTS, UnsafeWritePath, guard_write_path


class PathGuardTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_every_protected_root_and_descendant_is_rejected(self) -> None:
        for relative in PROTECTED_RELATIVE_ROOTS:
            with self.subTest(relative=relative), self.assertRaises(UnsafeWritePath):
                guard_write_path(self.root / relative / "would_write.bin", self.root)

    def test_every_v02_root_is_allowed(self) -> None:
        for relative in ALLOWED_RELATIVE_ROOTS:
            expected = (self.root / relative / "new.json").resolve()
            with self.subTest(relative=relative):
                self.assertEqual(guard_write_path(expected, self.root), expected)

    def test_prefix_collision_and_escape_are_rejected(self) -> None:
        with self.assertRaises(UnsafeWritePath):
            guard_write_path(self.root / "project/results/pacer_fkg_v02_longmd_v01_evil/x", self.root)
        with self.assertRaises(UnsafeWritePath):
            guard_write_path(self.root / "project/pacer_fkg_v02/../../../outside/x", self.root)

    def test_normalized_v02_descendant_is_allowed(self) -> None:
        value = self.root / "project/pacer_fkg_v02/a/../b/report.json"
        self.assertEqual(guard_write_path(value, self.root), (self.root / "project/pacer_fkg_v02/b/report.json").resolve())


if __name__ == "__main__":
    unittest.main()
