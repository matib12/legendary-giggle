from __future__ import annotations

import asyncio
import csv
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
COMPOSE_FILE = PROJECT_ROOT / "tests" / "integration" / "docker-compose.yml"
SRC_DIR = PROJECT_ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from influx_client import InfluxClientFactory
from query_builder import QueryBuilder


logger = logging.getLogger("tests.integration.influx")


def _docker_compose_command() -> list[str] | None:
    if shutil.which("docker") is None:
        return None
    return ["docker", "compose"]


def _docker_daemon_available() -> bool:
    if shutil.which("docker") is None:
        return False

    result = subprocess.run(
        ["docker", "info"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def _wait_for_influx_startup(config_path: Path, timeout_s: int = 60) -> str:
    deadline = time.time() + timeout_s
    last_error = ""
    logger.info("Waiting for InfluxDB startup (timeout=%ss)", timeout_s)

    while time.time() < deadline:
        try:
            client = InfluxClientFactory(str(config_path)).build()
            try:
                if client.ping():
                    version = str(client.version())
                    logger.info("InfluxDB ping ok, version=%s", version)
                    return version
            finally:
                client.close()
        except Exception as exc:
            last_error = str(exc)
            logger.debug("InfluxDB not ready yet: %s", last_error)

        time.sleep(2)

    raise RuntimeError(f"InfluxDB did not become ready in {timeout_s}s: {last_error}")


@pytest.fixture(scope="module")
def test_influxdb() -> None:
    if not _docker_daemon_available():
        pytest.skip("docker daemon is not running")

    compose_cmd = _docker_compose_command()
    if compose_cmd is None:
        pytest.skip("docker is not available")

    up_command = compose_cmd + ["-f", str(COMPOSE_FILE), "up", "-d"]
    down_command = compose_cmd + ["-f", str(COMPOSE_FILE), "down", "-v"]
    wait_seed_command = compose_cmd + [
        "-f",
        str(COMPOSE_FILE),
        "wait",
        "influxdb-seed",
    ]

    logger.info("Starting docker compose stack for integration tests")
    logger.debug("Compose file: %s", COMPOSE_FILE)
    subprocess.run(up_command, check=True, cwd=PROJECT_ROOT)
    logger.info("Waiting for dummy-data seed container completion")
    subprocess.run(wait_seed_command, check=True, cwd=PROJECT_ROOT)
    time.sleep(2)
    logger.info("Integration test stack is ready")
    yield
    logger.info("Stopping docker compose stack")
    subprocess.run(down_command, check=False, cwd=PROJECT_ROOT)


def test_end_to_end_query_and_csv(tmp_path: Path, test_influxdb: None) -> None:
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "influx_url": "http://localhost:18086",
                "bucket": "test-bucket",
                "org": "test-org",
                "token": "test-token",
            }
        ),
        encoding="utf-8",
    )
    logger.info("Generated temporary integration config: %s", config_path)

    version = _wait_for_influx_startup(config_path)
    assert version

    params_path = tmp_path / "params.json"
    params_path.write_text(
        json.dumps(
            {
                "bucket": "test-bucket",
                "measurement": "cpu",
                "field": "usage",
                "start": "2024-01-01T00:00:00Z",
                "stop": "2030-01-01T00:00:00Z",
                "extra_filters": [["host", "seeded"]],
            }
        ),
        encoding="utf-8",
    )
    logger.info("Generated query params file: %s", params_path)

    output_csv = tmp_path / "read_cpu.csv"
    client = InfluxClientFactory(str(config_path)).build()
    try:
        current_dir = os.getcwd()
        os.chdir(PROJECT_ROOT)
        try:
            query_builder = QueryBuilder(client)
            generated_csv = asyncio.run(
                query_builder.run_query(
                    "queries/base_query.flux",
                    str(params_path),
                    str(output_csv),
                )
            )
            logger.info("Query execution completed, CSV path: %s", generated_csv)
        finally:
            os.chdir(current_dir)
    finally:
        client.close()

    generated_csv_path = Path(generated_csv)
    assert generated_csv_path.exists()

    with generated_csv_path.open("r", encoding="utf-8", newline="") as csv_file:
        rows = list(csv.DictReader(csv_file))
    logger.info("CSV rows read: %d", len(rows))
    assert len(rows) > 0