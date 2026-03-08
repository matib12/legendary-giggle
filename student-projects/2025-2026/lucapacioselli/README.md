# influx-scraper

## Problem description

This project automates data extraction from InfluxDB using parameterized Flux templates.
The goal is to make query execution reproducible and testable by:

- loading connection settings from `config.json`;
- validating query parameters before rendering the template;
- executing one or more queries asynchronously;
- saving Influx responses to CSV files for downstream analyses (To Be Implemented);
- testing the full flow against a real InfluxDB instance started with Docker.

## Project structure

```text
lucapacioselli/
├── LICENSE
├── README.md
├── config.json.template
├── data/
│   └── csv/
├── queries/
│   ├── base_query.flux
│   └── params.json.template
├── requirements.txt
├── src/
│   ├── influx_client.py
│   ├── main.py
│   ├── query_builder.py
│   └── validate_schema.py
└── tests/
    ├── integration/
    │   ├── docker-compose.yml
    │   ├── dummy_data.lp
    │   └── test_influx_integration.py
    └── test_client.py
```

## Reproducible setup

### 1) Create and activate a virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 2) Install exact dependencies

```bash
pip install -r requirements.txt
```

### 3) Create local config from template

```bash
cp config.json.template config.json
```

Then edit `config.json` with your real Influx values.

## Code walkthrough

### `src/influx_client.py`

- `InfluxClientFactory` loads JSON config once and builds an authenticated `InfluxDBClient`.
- Centralizes connection parameters so the rest of the app does not duplicate setup logic.

### `src/validate_schema.py`

- Defines `FluxQueryParams` (Pydantic) to validate query inputs.
- Enforces identifier/time formats before template rendering.

### `src/query_builder.py`

- Loads params JSON, validates it, renders `queries/base_query.flux` with Jinja2.
- Executes Flux queries via Influx client and writes raw CSV output.
- Supports async orchestration (`run_query`, `run_multiple`).

### `src/main.py`

- CLI entrypoint for ping mode and query execution mode.
- Builds jobs for each params file and writes timestamped CSV outputs.

## Run the main code

### Ping InfluxDB

```bash
python src/main.py --config config.json --ping
```

### Execute multiple query parameter files (adapt params as needed)

```bash
python src/main.py \
	--config config.json \
	--template queries/base_query.flux \
	--params queries/params.json.template queries/params2.json.template \
	--output-dir data/csv \
	--log-level INFO
```

## Tests

### Unit tests

- `tests/test_client.py`: factory/config behavior and client creation error propagation.

Run:

```bash
python -m pytest -q tests/test_client.py
```

### Integration test (real InfluxDB with Docker)

- `tests/integration/docker-compose.yml` starts:
	- `influxdb` (v2.7, healthcheck);
	- `influxdb-seed` (loads `dummy_data.lp` automatically on startup).
- `tests/integration/test_influx_integration.py` verifies:
	1. InfluxDB ping/version,
	2. real query execution via `QueryBuilder`,
	3. generated CSV exists and contains rows.

Run:

```bash
python -m pytest -q tests/integration/test_influx_integration.py
```

Run with logs:

```bash
python -m pytest tests/integration/test_influx_integration.py -s -vv --log-cli-level=INFO
```

## Use Docker test stack with `main.py`

You can reuse the integration Docker stack to manually test the main code on seeded dummy data.

### 1) Start Influx + seed data

```bash
docker compose -f tests/integration/docker-compose.yml up -d
docker compose -f tests/integration/docker-compose.yml wait influxdb-seed
```

### 2) Create a temporary test config

```bash
cat > /tmp/config.test.json <<'JSON'
{
	"influx_url": "http://localhost:18086",
	"bucket": "test-bucket",
	"org": "test-org",
	"token": "test-token"
}
JSON
```

### 3) Run main against dummy data (adapt params as needed to match `dummy_data.lp`)

```bash
python src/main.py \
	--config /tmp/config.test.json \
	--template queries/base_query.flux \
	--params queries/params.json.template \
	--output-dir data/csv
```

### 4) Stop containers

```bash
docker compose -f tests/integration/docker-compose.yml down -v
```

## Notes

- `config.json.template` is the shareable template for collaborators.
- `queries/params.json.template` is a template for query parameters, which can be copied and edited for different queries.
- Some JSON params are edited in place by CLI when `--start`/`--stop` are passed.

## External resources / attribution

- InfluxDB Python client: https://github.com/influxdata/influxdb-client-python
- InfluxDB Docker image: https://hub.docker.com/_/influxdb
- Jinja2 templating engine: https://jinja.palletsprojects.com/
- Pydantic validation: https://docs.pydantic.dev/
- Pytest framework: https://docs.pytest.org/
