# ARCTERX Slocum Data Harvester

Grab Slocum glider data, parse it, collate, and produce CF-1.13 compliant NetCDF files
for ARCTERX 2023.

## Installation

Install as a Python package (provides `log-harvest`, `decode-argos`, and `mk-combined` commands):

```bash
pip install .
```

Or install globally via pipx:

```bash
pipx install .
```

## Required Tools

### External binaries

| Tool | Purpose | Source |
|------|---------|--------|
| `rsync` | Sync glider data from remote servers | System package manager |
| `dbd2netCDF` | Convert Slocum binary (.sbd/.tbd) files to NetCDF | git@github.com:mousebrains/dbd2netcdf.git |

### Python dependencies

Installed automatically via `pip install .`:

numpy, pandas, xarray, scipy, gsw, netcdf4, TPWUtils

## Commands

### syncit (shell script)

Wrapper that invokes `syncit.py` with default directories:

```bash
./syncit
```

Runs the full pipeline with `--tgt=gliders --output=data --combo=googleDrive`.

### syncit.py

Main orchestrator. Rsyncs glider data from remote servers, then runs each
processing stage in sequence. Not installed as a command — run directly from
the repository.

```
syncit.py [--hostname HOST] [--src SRC] [--glider GLIDER]
          [--t0 T0] [--tgt DIR] [--output DIR] [--combo DIR]
          [--rsync PATH] [--dbd2netCDF PATH]
          [--verbose | --debug]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--hostname` | `fs2` | Remote server to rsync from |
| `--src` | auto | Source paths for rsync (repeatable) |
| `--glider` | `osu684 osu685 osu686` | Glider IDs to process (repeatable) |
| `--t0` | `20230527T000000` | Earliest log file timestamp to consider |
| `--tgt` | `.` | Target directory for rsynced raw data |
| `--output` | `.` | Output directory for intermediate NetCDF files |
| `--combo` | `.` | Output directory for final combined NetCDF files |
| `--rsync` | `/usr/bin/rsync` | Path to rsync binary |
| `--dbd2netCDF` | `/usr/local/bin/dbd2netCDF` | Path to dbd2netCDF binary |

**Pipeline stages:**
1. Rsync glider data and ARGOS messages from remote servers
2. Decode ARGOS satellite messages → `argos.nc`
3. Harvest glider log files → `log.nc`
4. Convert flight binary files (`.s?d`) → `flt.{glider}.nc`
5. Convert science binary files (`.t?d`) → `sci.{glider}.nc`
6. Combine all sources → `{glider}.nc`

### log-harvest

Parses Slocum glider log files and extracts GPS positions, sensor readings,
and water velocity into a single NetCDF file.

```
log-harvest [--t0 T0] [--nc OUTPUT] [--verbose | --debug] [FILE ...]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--t0` | `20230527T000000` | Earliest timestamp to include |
| `--nc` | `log.nc` | Output NetCDF filename |
| `FILE` | `osu68?/logs/*.log` | Log files to parse (glob if omitted) |

### decode-argos

Decodes ARGOS satellite position messages into a NetCDF file.

```
decode-argos [--nc OUTPUT] FILE [FILE ...]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--nc` | `tpw.nc` | Output NetCDF filename |
| `FILE` | (required) | ARGOS message files to decode |

### mk-combined

Merges log, flight, and science data for a single glider into a CF-1.8
compliant NetCDF file with derived oceanographic variables (salinity,
potential temperature, density).

```
mk-combined --glider ID --output FILE
            [--prefix PREFIX] [--ncLog FILE]
            [--ncFlight FILE] [--ncScience FILE]
            [--verbose | --debug]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--glider` | (required) | Glider number (e.g., `684`) |
| `--output` | (required) | Output NetCDF filename |
| `--prefix` | `osu` | Institution prefix for glider name |
| `--ncLog` | `log.nc` | Input log NetCDF from log-harvest |
| `--ncFlight` | auto | Input flight NetCDF (default: `flt.{prefix}{glider}.nc`) |
| `--ncScience` | auto | Input science NetCDF (default: `sci.{prefix}{glider}.nc`) |

## Quick Start

```bash
pip install .
./syncit
```

Or run the pipeline manually:

```bash
./syncit.py --verbose --tgt=gliders --output=data --combo=googleDrive
```
