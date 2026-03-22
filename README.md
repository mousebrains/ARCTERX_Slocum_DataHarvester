# ARCTERX Slocum Data Harvester

Grab Slocum glider data, parse it, collate, and produce CF-1.8 compliant NetCDF files
for the ARCTERX 2023 In-situ Observation Project.

## Required Tools

### External binaries

| Tool | Purpose | Source |
|------|---------|--------|
| `rsync` | Sync glider data from remote servers | System package manager |
| `dbd2netCDF` | Convert Slocum binary (.sbd/.tbd) files to NetCDF | git@github.com:mousebrains/dbd2netcdf.git |

### Python packages

Install with `pip install -r requirements.txt`:

- **numpy** — numerical computations
- **pandas** — tabular data handling
- **xarray** — N-dimensional datasets and NetCDF I/O
- **scipy** — interpolation of GPS positions to science times
- **gsw** — Gibbs Sea Water oceanographic calculations (salinity, density, etc.)
- **netcdf4** — NetCDF file backend for xarray

### Vendored library

**TPWUtils** is a snapshot of git@github.com:mousebrains/TPWUtils.git
(logging, threading, credential management utilities).

## Commands

### syncit (shell script)

Wrapper that invokes `syncit.py` with default directories:

```bash
./syncit
```

Runs the full pipeline with `--tgt=gliders --output=data --combo=googleDrive`.

### syncit.py

Main orchestrator. Rsyncs glider data from remote servers, then runs each
processing stage in sequence.

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

### logHarvest.py

Parses Slocum glider log files and extracts GPS positions, sensor readings,
and water velocity into a single NetCDF file.

```
logHarvest.py [--t0 T0] [--nc OUTPUT] [--verbose | --debug] [FILE ...]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--t0` | `20230527T000000` | Earliest timestamp to include |
| `--nc` | `log.nc` | Output NetCDF filename |
| `FILE` | `osu68?/logs/*.log` | Log files to parse (glob if omitted) |

### decodeArgos.py

Decodes ARGOS satellite position messages into a NetCDF file.

```
decodeArgos.py [--nc OUTPUT] FILE [FILE ...]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--nc` | `tpw.nc` | Output NetCDF filename |
| `FILE` | (required) | ARGOS message files to decode |

### mkCombined.py

Merges log, flight, and science data for a single glider into a CF-1.8
compliant NetCDF file with derived oceanographic variables (salinity,
potential temperature, density).

```
mkCombined.py --glider ID --output FILE
              [--prefix PREFIX] [--ncLog FILE]
              [--ncFlight FILE] [--ncScience FILE]
              [--verbose | --debug]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--glider` | (required) | Glider number (e.g., `684`) |
| `--output` | (required) | Output NetCDF filename |
| `--prefix` | `osu` | Institution prefix for glider name |
| `--ncLog` | `log.nc` | Input log NetCDF from logHarvest |
| `--ncFlight` | auto | Input flight NetCDF (default: `flt.{prefix}{glider}.nc`) |
| `--ncScience` | auto | Input science NetCDF (default: `sci.{prefix}{glider}.nc`) |

### slocum_utils.py

Shared utility module (not a standalone command). Provides `mkDegrees_scalar()`
and `mkDegrees()` for converting Slocum DDMM.MM coordinates to decimal degrees.

## Quick Start

```bash
pip install -r requirements.txt
./syncit
```

Or run the pipeline manually:

```bash
./syncit.py --verbose --tgt=gliders --output=data --combo=googleDrive
```
