#!/usr/bin/env python3
"""Generate native WHAM input files from a compact JSON specification."""

import argparse
import json
import sys
from pathlib import Path


TIMESERIES_KEYS = [
    "path",
    "columns",
    "skipfrombeginning",
    "skip",
    "outputs",
    "outputNames",
    "verbose",
]

BIAS_KEYS = [
    "type",
    "dimension",
    "xstar",
    "kappa",
    "phi",
    "temperature",
]

WHAM_KEYS = [
    "name",
    "type",
    "strategy",
    "precision",
    "verbose",
]

STRATEGY_KEYS = [
    "type",
    "name",
    "max_iterations",
    "epsilon",
    "epsilon_rel",
    "tolerance",
    "printevery",
]

BIN_KEYS = [
    "dimension",
    "range",
    "numbins",
]


class SpecError(ValueError):
    pass


def format_scalar(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        raise SpecError("Cannot write null values to WHAM input")
    return str(value)


def format_value(value):
    if isinstance(value, list):
        return "[ " + " ".join(format_scalar(item) for item in value) + " ]"
    return format_scalar(value)


def write_key_value(lines, key, value, indent):
    lines.append(f"{' ' * indent}{key} = {format_value(value)}")


def write_block(lines, name, entries, indent=0):
    lines.append(f"{' ' * indent}{name} = {{")
    for key, value in entries:
        if isinstance(value, list) and value and all(isinstance(item, tuple) for item in value):
            write_block(lines, key, value, indent + 4)
        else:
            write_key_value(lines, key, value, indent + 4)
    lines.append(f"{' ' * indent}}}")


def ordered_entries(mapping, keys):
    return [(key, mapping[key]) for key in keys if key in mapping]


def require_mapping(spec, key):
    value = spec.get(key)
    if not isinstance(value, dict):
        raise SpecError(f"'{key}' must be an object")
    return value


def require_list(mapping, key):
    value = mapping.get(key)
    if not isinstance(value, list) or not value:
        raise SpecError(f"'{key}' must be a non-empty list")
    return value


def infer_dimension(bias, data):
    if "dimension" in bias:
        return bias["dimension"]

    for key in ("xstar", "kappa", "phi"):
        value = bias.get(key)
        if isinstance(value, list) and value:
            return len(value)

    columns = data.get("columns")
    if isinstance(columns, list) and columns:
        return len(columns)

    raise SpecError("Bias dimension could not be inferred")


def merge_window_timeseries(data, window):
    fields = {}
    for key in TIMESERIES_KEYS:
        if key in data:
            fields[key] = data[key]

    if isinstance(window.get("timeseries"), dict):
        fields.update(window["timeseries"])

    for key in TIMESERIES_KEYS:
        if key in window:
            fields[key] = window[key]

    if "path" not in fields:
        raise SpecError("Each window must provide a timeseries path")
    if "columns" not in fields:
        raise SpecError("Each window must provide columns, or data.columns must be set")

    return fields


def merge_window_bias(data, window):
    fields = {}
    bias_defaults = data.get("bias_defaults", {})
    if bias_defaults:
        if not isinstance(bias_defaults, dict):
            raise SpecError("'data.bias_defaults' must be an object")
        fields.update(bias_defaults)

    if isinstance(window.get("bias"), dict):
        fields.update(window["bias"])

    for key in BIAS_KEYS:
        if key in window:
            fields[key] = window[key]

    fields["dimension"] = infer_dimension(fields, data)

    if not any(key in fields for key in ("xstar", "kappa", "phi")):
        raise SpecError("Each window must provide bias parameters")

    return fields


def build_timeseries_and_bias_blocks(spec):
    data = require_mapping(spec, "data")
    windows = require_list(data, "windows")

    blocks = []
    for window in windows:
        if not isinstance(window, dict):
            raise SpecError("Each data window must be an object")
        blocks.append(("timeseries", ordered_entries(merge_window_timeseries(data, window), TIMESERIES_KEYS)))

    for window in windows:
        blocks.append(("bias", ordered_entries(merge_window_bias(data, window), BIAS_KEYS)))

    return blocks


def build_outputs(wham):
    outputs = wham.get("outputs")
    if outputs is None:
        return []

    if isinstance(outputs, dict):
        return [
            ("outputs", list(outputs.keys())),
            ("outputFile", list(outputs.values())),
        ]

    if isinstance(outputs, list):
        output_files = wham.get("outputFile")
        if not isinstance(output_files, list):
            raise SpecError("'wham.outputFile' must be provided when 'wham.outputs' is a list")
        if len(outputs) != len(output_files):
            raise SpecError("'wham.outputs' and 'wham.outputFile' must have the same length")
        return [
            ("outputs", outputs),
            ("outputFile", output_files),
        ]

    raise SpecError("'wham.outputs' must be an object or a list")


def build_strategy_entries(wham):
    strategies = wham.get("strategies")
    if strategies is None:
        return []

    if not isinstance(strategies, list) or not strategies:
        raise SpecError("'wham.strategies' must be a non-empty list")

    entries = []
    strategy_names = []
    for strategy in strategies:
        if not isinstance(strategy, dict):
            raise SpecError("Each WHAM strategy must be an object")
        if "type" not in strategy:
            raise SpecError("Each WHAM strategy must include a type")
        if "name" not in strategy:
            raise SpecError("Each WHAM strategy must include a name")
        entries.append(("Uwhamstrategy", ordered_entries(strategy, STRATEGY_KEYS)))
        strategy_names.append(strategy["name"])

    explicit_names = wham.get("strategyNames")
    if explicit_names is not None:
        if not isinstance(explicit_names, list) or not explicit_names:
            raise SpecError("'wham.strategyNames' must be a non-empty list")
        strategy_names = explicit_names

    entries.append(("strategyNames", strategy_names))

    return entries


def build_bin_entries(wham):
    bins = require_list(wham, "bins")
    entries = []
    for bin_spec in bins:
        if not isinstance(bin_spec, dict):
            raise SpecError("Each bin specification must be an object")
        if "range" not in bin_spec or "numbins" not in bin_spec:
            raise SpecError("Each bin specification must include range and numbins")
        entries.append(("bins", ordered_entries(bin_spec, BIN_KEYS)))
    return entries


def build_wham_block(spec):
    wham = require_mapping(spec, "wham")
    if "type" not in wham:
        raise SpecError("'wham.type' is required")

    entries = ordered_entries(wham, WHAM_KEYS)

    if wham["type"] == "Uwham":
        if "strategies" not in wham:
            raise SpecError("'wham.strategies' is required when wham.type is Uwham")
        entries.extend(build_strategy_entries(wham))

    entries.extend(build_bin_entries(wham))

    for key in ("ErrorAnalysis", "ErrorIteration", "BAR", "N", "Nvec"):
        if key in wham:
            entries.append((key, wham[key]))

    entries.extend(build_outputs(wham))

    return ("wham", entries)


def generate(spec):
    blocks = build_timeseries_and_bias_blocks(spec)
    blocks.append(build_wham_block(spec))

    lines = []
    for index, (name, entries) in enumerate(blocks):
        if index:
            lines.append("")
        write_block(lines, name, entries)

    return "\n".join(lines) + "\n"


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spec", type=Path, help="JSON WHAM input specification")
    parser.add_argument("-o", "--output", type=Path, help="Output path. Defaults to stdout.")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        with args.spec.open() as handle:
            spec = json.load(handle)
        output = generate(spec)
    except (OSError, json.JSONDecodeError, SpecError) as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 1

    if args.output:
        args.output.write_text(output)
    else:
        sys.stdout.write(output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
