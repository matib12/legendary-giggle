from __future__ import annotations

import json
from pathlib import Path

import pytest

from src.influx_client import InfluxClientFactory


def _write_config(path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "influx_url": "http://localhost:8086",
                "bucket": "test-bucket",
                "org": "test-org",
                "token": "test-token",
            }
        ),
        encoding="utf-8",
    )


def test_factory_loads_config_from_json(tmp_path: Path) -> None:
    config_path = tmp_path / "config.json"
    _write_config(config_path)

    factory = InfluxClientFactory(str(config_path))

    assert factory.config["influx_url"] == "http://localhost:8086"
    assert factory.config["bucket"] == "test-bucket"


def test_factory_raises_missing_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing.json"

    with pytest.raises(FileNotFoundError):
        InfluxClientFactory(str(missing))


def test_factory_build_creates_influx_client(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "config.json"
    _write_config(config_path)

    captured: dict[str, str] = {}

    class FakeInfluxDBClient:
        def __init__(self, *, url: str, token: str, org: str):
            captured["url"] = url
            captured["token"] = token
            captured["org"] = org

    monkeypatch.setattr("src.influx_client.InfluxDBClient", FakeInfluxDBClient)

    factory = InfluxClientFactory(str(config_path))
    client = factory.build()

    assert isinstance(client, FakeInfluxDBClient)
    assert captured == {
        "url": "http://localhost:8086",
        "token": "test-token",
        "org": "test-org",
    }


def test_factory_build_propagates_library_errors(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "config.json"
    _write_config(config_path)

    class BoomError(RuntimeError):
        pass

    def _raise_error(**_: str):
        raise BoomError("cannot connect")

    monkeypatch.setattr("src.influx_client.InfluxDBClient", _raise_error)

    factory = InfluxClientFactory(str(config_path))

    with pytest.raises(BoomError):
        factory.build()