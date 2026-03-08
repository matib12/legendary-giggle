from typing import ClassVar, Optional, List
from pydantic import BaseModel, field_validator, ValidationError
import re


class FluxQueryParams(BaseModel):
    """Validated input model used before rendering Flux queries."""

    bucket: str
    measurement: str
    field: str
    start: str
    stop: str
    extra_filters: Optional[List[List[str]]] = None
    aggregate: Optional[str] = None

    IDENTIFIER_REGEX: ClassVar[re.Pattern] = re.compile(r"^[a-zA-Z0-9_\-]+$")
    TIME_REGEX: ClassVar[re.Pattern] = re.compile(
        r"^(-?\d+[smhdw]|now\(\)|\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(\.\d+)?Z)$"
    )

    @field_validator("bucket", "measurement", "field")
    @classmethod
    def validate_identifier(cls, v: str) -> str:
        # Restrict identifiers to safe alphanumeric/tag-friendly characters.
        if not cls.IDENTIFIER_REGEX.match(v):
            raise ValueError(f"Invalid identifier: {v}")
        return v

    @field_validator("start", "stop")
    @classmethod
    def validate_time(cls, v: str) -> str:
        # Accept relative durations (e.g. -1h), now(), or full RFC3339 timestamps.
        if not cls.TIME_REGEX.match(v):
            raise ValueError(f"Invalid time format: {v}")
        return v

    @field_validator("extra_filters")
    @classmethod
    def validate_filters(cls, v: Optional[List[List[str]]]) -> Optional[List[List[str]]]:
        # Require filters as [field, value] pairs.
        if v is None:
            return v

        for item in v:
            if not isinstance(item, list) or len(item) != 2:
                raise ValueError(f"Invalid filter: {item}")
        return v


def validate_params(data: dict) -> FluxQueryParams:
    """Validate raw JSON params and return a typed model."""
    return FluxQueryParams(**data)
    