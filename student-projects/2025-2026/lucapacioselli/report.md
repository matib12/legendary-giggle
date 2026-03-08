# Performance Report

## Scope

This report summarizes one **example debug run** of the main pipeline with two query parameter files:

```json
'cpu_usage.json': {
    "bucket": "test-bucket",
    "start": "2024-01-01T00:00:00Z",
    "stop": "2030-01-01T00:00:00Z",
    "measurement": "cpu",
    "field": "usage",
    "extra_filters": [
        ["host", "seeded"]
    ],
    "aggregate": "aggregateWindow(every: 10s, fn: mean, createEmpty: false)"
 }
```

```json
 'mem_usage.json': {
    "bucket": "test-bucket",
    "start": "2024-01-01T00:00:00Z",
    "stop": "2030-01-01T00:00:00Z",
    "measurement": "memory",
    "field": "used_percent",
    "extra_filters": [
        ["host", "seeded"]
    ],
    "aggregate": "aggregateWindow(every: 10s, fn: mean, createEmpty: false)"
 }
 ```

The run was executed with `time.perf_counter()` instrumentation and `cProfile` enabled in `DEBUG` mode.

## Run setup

- Python: `.venv` interpreter
- Data source: test InfluxDB from `tests/integration/docker-compose.yml` (with seeded dummy data)
- Config file: `/tmp/lucapacioselli_profile_run/config.test.json`

```json
{
    "influx_url": "http://localhost:18086",
    "bucket": "test-bucket",
    "org": "test-org",
    "token": "test-token"
}
```

- Params files:
  - `/tmp/lucapacioselli_profile_run/cpu_usage.json`
  - `/tmp/lucapacioselli_profile_run/mem_usage.json`
- Command used:

```bash
python src/main.py \
  --config /tmp/lucapacioselli_profile_run/config.test.json \
  --template queries/base_query.flux \
  --params /tmp/lucapacioselli_profile_run/cpu_usage.json /tmp/lucapacioselli_profile_run/mem_usage.json \
  --output-dir /tmp/lucapacioselli_profile_run/out \
  --start 2024-01-01T00:00:00Z \
  --stop 2030-01-01T00:00:00Z \
  --log-level DEBUG
```

## Outputs generated

- `/tmp/lucapacioselli_profile_run/out/cpu_usage_2026-03-08_13-21-50.csv`
- `/tmp/lucapacioselli_profile_run/out/mem_usage_2026-03-08_13-21-50.csv`

## Timing results (`perf_counter`)

Captured by `main.py` in debug mode:

- `setup_client_s`: `0.006896`
- `prepare_jobs_s`: `0.000971`
- `run_queries_s`: `0.068649`
- `postprocess_logs_s`: `0.000120`
- `total_s`: `0.076731`

Per-query debug timings from `query_builder.py`:

- `cpu_usage.json`: `0.068430s`
- `mem_usage.json`: `0.061498s`

## Profiling results (`cProfile`)

Summary:

- `20352` function calls (`19329` primitive)
- total profiled time: `0.077s`

Top cumulative-time hotspots indicate the dominant cost is network I/O in InfluxDB HTTP calls:

1. `influxdb_client.client.query_api.query_csv`
2. `influxdb_client.service.query_service.post_query`
3. `urllib3` request/connection stack (`urlopen`, `_make_request`, `getresponse`)
4. socket/http read operations (`http.client`, `socket.readinto`, buffered `readline`)

## Bottleneck interpretation

For this sample run, the main bottleneck is **remote query execution and HTTP response read**, not local Python processing.

- Local overhead (job preparation + logging) is negligible.
- End-to-end latency is mostly attributable to InfluxDB round-trip and response streaming.

## Notes

- Results are from a local Docker-based test environment; absolute times may differ on production datasets/network.
