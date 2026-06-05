#!/usr/bin/env python3
"""Fail if any tracked text artifact contains CJK characters (English-only repo).

External/public artifacts must read as natural English. This scans tracked text
files for CJK ranges (CJK punctuation, Hiragana/Katakana, CJK ideographs, Hangul,
fullwidth forms) and exits non-zero if any are found. The character class is built
from integer code points so this checker's own source stays pure ASCII and never
flags itself.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

# (low, high) inclusive code-point ranges, built via chr() to keep this file ASCII.
_RANGES = [
    (0x3000, 0x303F),  # CJK symbols & punctuation
    (0x3040, 0x30FF),  # Hiragana + Katakana
    (0x3400, 0x4DBF),  # CJK extension A
    (0x4E00, 0x9FFF),  # CJK unified ideographs
    (0xAC00, 0xD7AF),  # Hangul syllables
    (0xFF00, 0xFFEF),  # fullwidth / halfwidth forms
]
CJK = re.compile("[" + "".join(f"{chr(lo)}-{chr(hi)}" for lo, hi in _RANGES) + "]")

TEXT_EXT = {".py", ".md", ".yaml", ".yml", ".toml", ".sh", ".json", ".txt", ".cfg", ".ini", ".rst"}


def _tracked_files() -> list[Path]:
    try:
        out = subprocess.run(
            ["git", "ls-files"], capture_output=True, text=True, check=True
        ).stdout
        files = [Path(p) for p in out.splitlines() if p]
        if files:
            return files
    except (subprocess.SubprocessError, OSError):
        pass
    return [p for p in Path().rglob("*") if p.is_file()]


def main() -> int:
    bad: list[str] = []
    scanned = 0
    for p in _tracked_files():
        if p.suffix.lower() not in TEXT_EXT or not p.exists():
            continue
        scanned += 1
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        if CJK.search(text):
            bad.append(str(p))

    if bad:
        print("english-only: FAIL -- CJK characters found in:")
        for b in bad:
            print(f"  {b}")
        return 1
    print(f"english-only: OK -- scanned {scanned} text artifacts; no CJK characters found.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
