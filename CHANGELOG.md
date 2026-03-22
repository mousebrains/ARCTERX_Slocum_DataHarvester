# Changelog

## [1.0.0] - 2026-03-22

### Added
- Installable Python package (`pip install .` / `pipx install .`)
- CLI entry points: `log-harvest`, `decode-argos`, `mk-combined`
- `pyproject.toml` with dependencies and console_scripts
- `requirements.txt`
- Pre-commit hooks (trailing whitespace, syntax checks)
- Comprehensive README with command documentation
- Input file and variable validation in `mkCombo`
- Empty dataset guards in `logHarvest` and `decodeArgos`
- Data reduction logging in `mkCombined`
- External tool validation (`rsync`, `dbd2netCDF`) at startup

### Changed
- CF compliance updated from 1.8 to 1.13
  - `depth` standard_name corrected to `depth` (was invalid `sea_water_depth`)
  - `rho` standard_name corrected to `sea_water_sigma_t` (stores density - 1000)
  - Added `units_metadata = "temperature: on_scale"` per CF-1.11
  - Added `featureType = "trajectory"` global attribute
  - `sigma` metadata clarified as potential density anomaly (sigma-0)
- `log-harvest` optimized: 8.5x faster via dict lookup and compiled regexes
- TPWUtils replaced with PyPI dependency (removed vendored copy)
- Modules moved into `slocum_harvester/` package
- `mkCombo` returns `bool` so callers can detect failure
- `syncit.py` main logic wrapped in `if __name__` guard
- Failed subprocess commands logged at WARNING with return code
- Metadata timestamps use UTC instead of local time
- Shell script (`syncit`) uses `set -euo pipefail`
- `interp1d` calls use `fill_value=np.nan` instead of string `"nan"`
- Shared `mkDegrees` extracted to `slocum_utils.py`; fixed vectorized sign handling

### Fixed
- Bare `except:` clauses replaced with `except Exception:`
- ARGOS timestamps now include UTC timezone
- `SettingWithCopyWarning` in `mkCombined` fixed with `.loc[]`
- `not currTime` check replaced with `currTime is None` (epoch-safe)
- Redundant `not currTime` removed from sensor parsing
- `interp1d` guarded against fewer than 2 data points
- Empty log.nc no longer crashes `mkCombo` (validates required variables)
- NaN lat/lon dropped before GSW oceanographic calculations
- Pipeline exits non-zero on unhandled exceptions
- Malformed datetime strings caught in log and ARGOS parsing

## [0.1.0] - 2023-06-29 to 2023-11-07

Initial development for ARCTERX 2023.

- `syncit.py` orchestrator with rsync, dbd2netCDF, and data combination
- `logHarvest.py` for parsing Slocum glider log files
- `decodeArgos.py` for ARGOS satellite message decoding
- `mkCombined.py` for producing combined CF-1.8 NetCDF output
- Vendored TPWUtils library
