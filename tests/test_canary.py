"""Canary + audit-chain smoke tests (substrate sanity, no torch/data needed)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mocnv import audit, canary  # noqa: E402


def test_canary_ok():
    result = canary.check()
    assert result["ok"] is True
    assert result["job_id"] == "canary-mocnv-demo"


def test_audit_chain_verifies(tmp_path):
    led = tmp_path / "led.ndjson"
    audit.emit("a", "job", {"x": 1}, ledger_path=led)
    audit.emit("b", "job", {"y": 2}, ledger_path=led)
    ok, n, bad = audit.verify(led)
    assert ok and n == 2 and bad is None


def test_audit_detects_tampering(tmp_path):
    led = tmp_path / "led.ndjson"
    audit.emit("a", "job", {"x": 1}, ledger_path=led)
    audit.emit("b", "job", {"y": 2}, ledger_path=led)
    lines = led.read_text().splitlines()
    entry = json.loads(lines[0])
    entry["fields"]["x"] = 999
    lines[0] = json.dumps(entry, sort_keys=True, separators=(",", ":"))
    led.write_text("\n".join(lines) + "\n")
    ok, _, _ = audit.verify(led)
    assert ok is False
