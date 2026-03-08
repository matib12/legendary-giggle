import json
import logging
from influxdb_client import InfluxDBClient

logger = logging.getLogger("InfluxClient")


class InfluxClientFactory:
    """Factory responsible for building an authenticated InfluxDB client.

    The class reads all connection parameters from `config.json` (or another
    JSON path passed by CLI) and keeps the parsed configuration in memory.
    """

    def __init__(self, config_path: str):
        # Load and cache configuration once at startup.
        self.config = self._load_json(config_path)

    def _load_json(self, path):
        # Any invalid/missing file errors are intentionally propagated.
        with open(path) as f:
            return json.load(f)

    def build(self) -> InfluxDBClient:
        """Create a new InfluxDB client using config values."""
        influx_cfg = self.config

        return InfluxDBClient(
            url=influx_cfg["influx_url"],
            token=influx_cfg["token"],
            org=influx_cfg["org"]
        )
