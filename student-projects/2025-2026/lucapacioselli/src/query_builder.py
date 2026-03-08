import asyncio
import json
import logging
import time
from jinja2 import Environment, FileSystemLoader, StrictUndefined
from validate_schema import validate_params, ValidationError

logger = logging.getLogger("QueryBuilder")


class QueryBuilder:
    """Render Flux templates and execute Influx queries asynchronously.

    Query execution runs in a worker thread (`asyncio.to_thread`) so multiple
    jobs can be scheduled concurrently from the async CLI.
    """

    def __init__(self, influx_client):
        self.client = influx_client
        self.query_api = self.client.query_api()

        self.jinja_env = Environment(
            loader=FileSystemLoader("."),
            undefined=StrictUndefined,
            trim_blocks=True,
            lstrip_blocks=True,
        )

    def _load_json(self, path):
        """Load query parameters from a JSON file."""
        with open(path) as f:
            logger.debug("Loading JSON from %s", path)
            return json.load(f)

    def _render_query(self, template_name: str, params: dict) -> str:
        """Validate parameters and render the Flux template."""
        validated = validate_params(params)

        template = self.jinja_env.get_template(template_name)
        rendered = template.render(**validated.model_dump())

        logger.debug("Rendered query:\n%s", rendered)
        return rendered

    def _execute_query(self, flux_query: str, output_path: str):
        """Execute a Flux query and persist raw CSV rows to disk."""
        csv_result = self.query_api.query_csv(flux_query)

        with open(output_path, "w") as f:
            for row in csv_result:
                f.write(",".join(row) + "\n")

        return output_path

    async def run_query(self, template_name: str, params_path: str, output_path: str):
        """Run one query job end-to-end and return the generated CSV path."""
        query_start = time.perf_counter()
        try:
            params = self._load_json(params_path)
            flux_query = self._render_query(template_name, params)

            result = await asyncio.to_thread(
                self._execute_query,
                flux_query,
                output_path
            )
            if logger.isEnabledFor(logging.DEBUG):
                elapsed = time.perf_counter() - query_start
                logger.debug("Query completed for %s in %.6fs", params_path, elapsed)
            return result
        except ValidationError as e:
            logger.error("Validation error: %s", e)
            raise

    async def run_multiple(self, jobs: list):
        """Schedule and await multiple query jobs concurrently."""
        tasks = [
            self.run_query(job["template"], job["params"], job["output"])
            for job in jobs
        ]
        return await asyncio.gather(*tasks)