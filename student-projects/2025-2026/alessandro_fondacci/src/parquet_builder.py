#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
C-V txt folder -> 1 row per file -> Parquet

Parquet columns:
Produzione, Wafer1, Wafer2, Shot1, Shot2, Struttura,
Temperatura, Frequenza, ACampl, Data, Ora, V, C

Filename example:
NLGAD1_W4_NtD_4-4_SP1_LGAD_20C_100Hz_100mV_24-02-2026_10h3m34s.txt

Wafer2 allowed: NtD | Deep | Shallow+Deep
Shot2 allowed:  SP1 | SP2

V = 1st column inside txt
C = 2nd column inside txt

Dependencies:
  pip install pandas pyarrow
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Tuple, Optional

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


FILENAME_RE = re.compile(
    r"^(?P<Produzione>[^_]+)_"
    r"W(?P<Wafer1>\d+)_"
    r"(?P<Wafer2>NtD|Deep|Shallow\+Deep)_"
    r"(?P<Shot1>[^_]+)_"
    r"(?P<Shot2>SP1|SP2)_"
    r"(?P<Struttura>[^_]+)_"
    r"(?P<TempVal>[-+]?\d+(?:\.\d+)?)C_"
    r"(?P<FreqVal>\d+(?:\.\d+)?)(?P<FreqUnit>Hz|kHz|MHz)_"
    r"(?P<ACVal>\d+(?:\.\d+)?)(?P<ACUnit>mV|V)_"
    r"(?P<d>\d{1,2})-(?P<m>\d{1,2})-(?P<y>\d{4})_"
    r"(?P<hh>\d{1,2})h(?P<mm>\d{1,2})m(?P<ss>\d{1,2})s\.txt$"
)


def convert_to_hz(val: float, unit: str) -> float:
    if unit == "Hz":
        return val
    if unit == "kHz":
        return val * 1e3
    if unit == "MHz":
        return val * 1e6
    raise ValueError(f"Unsupported frequency unit: {unit}")


def convert_to_volt(val: float, unit: str) -> float:
    if unit == "V":
        return val
    if unit == "mV":
        return val * 1e-3
    raise ValueError(f"Unsupported AC amplitude unit: {unit}")


def read_v_c_columns(txt_path: Path) -> Tuple[List[float], List[float]]:
    """
    Reads the first two columns from a tab-delimited txt file.
    Robust to extra columns and headers.
    """
    # Read only the first two columns; let pandas infer header.
    # - If there is a header row: pandas uses it and reads numeric below.
    # - If there isn't: header=0 will treat first row as header (bad). So we try header=0 first,
    #   and if we get no numeric data, retry with header=None.
    def _try_read(header):
        df = pd.read_csv(
            txt_path,
            sep="\t",
            engine="python",
            usecols=[0, 1],
            header=header,
            comment=None,
        )
        # Force numeric conversion (strings -> NaN)
        v = pd.to_numeric(df.iloc[:, 0], errors="coerce")
        c = pd.to_numeric(df.iloc[:, 1], errors="coerce")
        good = ~(v.isna() | c.isna())
        return v[good].astype(float).tolist(), c[good].astype(float).tolist()

    V, C = _try_read(header=0)
    if len(V) == 0:
        V, C = _try_read(header=None)

    if len(V) == 0:
        raise ValueError(f"No valid numeric (V,C) data found in: {txt_path}")

    return V, C


def txt_folder_to_parquet_cv(folder_path: str | Path, out_parquet_path: Optional[str | Path] = None) -> Path:
    folder = Path(folder_path).expanduser().resolve()
    if out_parquet_path is None or str(out_parquet_path).strip() == "":
        out_parquet = folder / "output_CV.parquet"
    else:
        out_parquet = Path(out_parquet_path).expanduser().resolve()

    txt_files = sorted(folder.glob("*.txt"))
    if not txt_files:
        raise FileNotFoundError(f"No .txt files found in: {folder}")

    rows = []
    for p in txt_files:
        m = FILENAME_RE.match(p.name)
        if not m:
            raise ValueError(f"Filename does not match expected pattern: {p.name}")

        g = m.groupdict()

        temperatura = float(g["TempVal"])
        frequenza_hz = convert_to_hz(float(g["FreqVal"]), g["FreqUnit"])
        acampl_v = convert_to_volt(float(g["ACVal"]), g["ACUnit"])

        data_str = f"{int(g['d']):02d}/{int(g['m']):02d}/{int(g['y']):04d}"
        ora_str = f"{int(g['hh']):02d}:{int(g['mm']):02d}:{int(g['ss']):02d}"

        V, C = read_v_c_columns(p)

        rows.append(
            dict(
                Produzione=g["Produzione"],
                Wafer1=int(g["Wafer1"]),
                Wafer2=g["Wafer2"],
                Shot1=g["Shot1"],
                Shot2=g["Shot2"],
                Struttura=g["Struttura"],
                Temperatura=temperatura,
                Frequenza=frequenza_hz,
                ACampl=acampl_v,
                Data=data_str,
                Ora=ora_str,
                V=V,  # list[float]
                C=C,  # list[float]
            )
        )

    df = pd.DataFrame(rows)

    # Write Parquet with explicit schema so V and C are true list columns (not generic objects)
    schema = pa.schema(
        [
            ("Produzione", pa.string()),
            ("Wafer1", pa.int64()),
            ("Wafer2", pa.string()),
            ("Shot1", pa.string()),
            ("Shot2", pa.string()),
            ("Struttura", pa.string()),
            ("Temperatura", pa.float64()),
            ("Frequenza", pa.float64()),
            ("ACampl", pa.float64()),
            ("Data", pa.string()),
            ("Ora", pa.string()),
            ("V", pa.list_(pa.float64())),
            ("C", pa.list_(pa.float64())),
        ]
    )

    table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)
    pq.write_table(table, out_parquet)

    print(f"Created Parquet: {out_parquet}")
    print(f"Rows written: {len(df)}")
    return out_parquet


if __name__ == "__main__":
    # Example usage:
    # python3 cvparquet.py "/Users/alefond/Desktop/NLGAD_CV_NotIrr"
    import sys

    if len(sys.argv) < 2:
        raise SystemExit("Usage: python3 cvparquet.py <folder_path> [out_parquet_path]")

    folder = sys.argv[1]
    outp = sys.argv[2] if len(sys.argv) >= 3 else None
    txt_folder_to_parquet_cv(folder, outp)