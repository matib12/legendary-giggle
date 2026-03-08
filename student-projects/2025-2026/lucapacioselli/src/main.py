import argparse
import asyncio
import cProfile
import json
import pstats
import time
import logging
from io import StringIO
from influx_client import InfluxClientFactory
from query_builder import QueryBuilder

logger = logging.getLogger("Main")


def ping_influx(client):
    """Ping InfluxDB and return server version when reachable."""
    try:
        if client.ping():
            version = client.version()
            return version
        return None
    except Exception as e:
        logger.error("Ping failed: %s", e)
        return None


def get_parser():
    """Build the CLI parser for ping mode and query execution mode."""
    parser = argparse.ArgumentParser(
        description="Async InfluxDB Query Runner"
    )

    parser.add_argument(
        "--config",
        default="config.json",
        help="Path to config.json"
    )

    parser.add_argument(
        "--template",
        default="queries/base_query.flux",
        help="Flux template name"
    )

    parser.add_argument(
        "--params",
        nargs="+",
        help="One or more params JSON files"
    )

    parser.add_argument(
        "--output-dir",
        default="data/csv",
        help="Directory where CSV files will be saved"
    )

    parser.add_argument(
        "--ping",
        action="store_true",
        help="Ping InfluxDB and print version"
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )

    parser.add_argument(
        "--start",
        default="-1h",
        help="Start time for the query (e.g., -30m, 2024-01-01T00:00:00Z). Default: -1h"
    )

    parser.add_argument(
        "--stop",
        default="now()",
        help="Stop time for the query (e.g., 2024-01-01T01:00:00Z). Default: now()"
    )

    return parser


async def main():
    """Run the CLI application.

    Flow:
    1) Build and validate Influx connection
    2) Optional ping mode
    3) Render and execute one or more Flux queries
    4) Save CSV outputs
    """

    parser = get_parser()
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    )

    debug_mode = logger.isEnabledFor(logging.DEBUG)
    profiler = cProfile.Profile() if debug_mode else None
    phase_timings: dict[str, float] = {}
    total_start = time.perf_counter()

    if profiler is not None:
        profiler.enable()

    setup_start = time.perf_counter()
    factory = InfluxClientFactory(args.config)
    influx_client = factory.build()
    phase_timings["setup_client_s"] = time.perf_counter() - setup_start

    if args.ping:
        ping_start = time.perf_counter()
        version = ping_influx(influx_client)
        phase_timings["ping_s"] = time.perf_counter() - ping_start
        print("InfluxDB version:", version)
        influx_client.close()
        total_elapsed = time.perf_counter() - total_start
        phase_timings["total_s"] = total_elapsed

        if debug_mode:
            print("\n=== Timing (seconds) ===")
            for name, elapsed in phase_timings.items():
                print(f"{name}: {elapsed:.6f}")

        if profiler is not None:
            profiler.disable()
            stream = StringIO()
            stats = pstats.Stats(profiler, stream=stream).sort_stats("cumulative")
            stats.print_stats(25)
            print("\n=== cProfile (top 25 by cumulative time) ===")
            print(stream.getvalue())
        return

    query_builder = QueryBuilder(influx_client)

    jobs = []
    jobs_start = time.perf_counter()

    for params_file in args.params:
        output_path = (
            f"{args.output_dir}/{params_file.split('/')[-1].replace('.json', '')}"
            + "_"
            + f"{time.strftime('%Y-%m-%d_%H-%M-%S')}.csv"
        )


        if args.start and args.stop:
            logger.info("Using start time: %s", args.start)
            logger.info("Using stop time: %s", args.stop)
            # Keep params files in sync with CLI time-window overrides.
            with open(params_file, "r") as f:
                params_dict = json.load(f)
                params_dict.update({
                    "start": args.start,
                    "stop": args.stop
                })
            with open(params_file, "w") as f:
                json.dump(params_dict, f, indent=4)

        jobs.append({
            "template": args.template,
            "params": params_file,
            "output": output_path
        })
    phase_timings["prepare_jobs_s"] = time.perf_counter() - jobs_start

    execute_start = time.perf_counter()
    results = await query_builder.run_multiple(jobs)
    phase_timings["run_queries_s"] = time.perf_counter() - execute_start

    save_log_start = time.perf_counter()
    for path in results:
        logger.info("Saved CSV to %s", path)
    phase_timings["postprocess_logs_s"] = time.perf_counter() - save_log_start

    influx_client.close()

    total_elapsed = time.perf_counter() - total_start
    phase_timings["total_s"] = total_elapsed

    if debug_mode:
        print("\n=== Timing (seconds) ===")
        for name, elapsed in phase_timings.items():
            print(f"{name}: {elapsed:.6f}")

    if profiler is not None:
        profiler.disable()
        stream = StringIO()
        stats = pstats.Stats(profiler, stream=stream).sort_stats("cumulative")
        stats.print_stats(25)
        print("\n=== cProfile (top 25 by cumulative time) ===")
        print(stream.getvalue())


if __name__ == "__main__":
    asyncio.run(main())
