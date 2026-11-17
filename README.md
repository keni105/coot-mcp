# coot-mcp

MCP servers that let an AI assistant (e.g. Claude Code) build and validate crystallographic models with
[Coot](https://github.com/pemsley/coot), both headless and in a running Coot window.
Tested with Coot 1.3.4 as shipped in Phenix 2.3rc2 on macOS.

| Server | What it does | Runs |
|---|---|---|
| `coot-edit` | `coot_edit_server.py` (this repo): headless model building and validation, 35 tools | Phenix's `phenix.python` |
| `coot-live` | drives a running Coot GUI through its JSON-RPC listener | Phenix's own `coot_mcp_socket_bridge.py`; start Coot with `coot-live` (this repo) |

Phenix also ships a read-only density-probe server (`phenix.coot_headless_mcp_server`) that works alongside these.

## Requirements

- Phenix 2.3 or later (it ships Coot 1.3.4 with `coot_headless_api`, and `fastmcp`).
- Set `PHENIX` to your Phenix install if it is not `/Applications/phenix-2.3rc2-6237`.

## Register with Claude Code

```bash
PHENIX=/Applications/phenix-2.3rc2-6237
claude mcp add --scope user coot-edit -- "$PHENIX/bin/phenix.python" "$PWD/coot_edit_server.py"
claude mcp add --scope user coot-live -- "$PHENIX/bin/python" -I "$PHENIX/share/coot/mcp/coot_mcp_socket_bridge.py"
```

## coot-edit

Load → inspect → edit → write:

- `coot_load_model`, `coot_load_maps` (phenix.refine labels by default), `coot_load_dictionary`, `coot_get_monomer`
- Validation: `coot_validate` (density fit/correlation, Q-score, rotamer, Rama, omega), `coot_geometry_outliers`,
  `coot_difference_peaks` (σ), `coot_unmodelled_blobs`, `coot_check_waters`, `coot_clashes`, `coot_residue_info`,
  `coot_residues_near`, `coot_density_at`
- Editing: `coot_refine` (SINGLE…ALL or LITERAL range), `coot_mutate`, `coot_auto_fit_rotamer`, `coot_flip_peptide`,
  `coot_side_chain_180`, `coot_add_alt_conf`, `coot_fill_partial_residues`, `coot_add_waters`, `coot_add_terminal_residue`,
  `coot_rigid_body_fit`, `coot_fit_ligand`, `coot_merge`, `coot_delete`, `coot_set_occupancy_b`, `coot_undo`, `coot_copy_model`
- `coot_write_model` (never overwrites unless `overwrite=True`)
- Escape hatch: `coot_api_search`, `coot_run_python` (`mc` = `molecules_container_t`)

Every edit reports the residues the CID actually matched and per-residue density (2mFo-DFc mean/min σ,
mFo-DFc min/max σ) before and after. Single-residue tools refuse CIDs that match more than one residue.

## coot-live

```bash
./coot-live                     # or: ./coot-live --pdb model.pdb --auto maps.mtz
```

Starts Phenix's Coot GUI with its JSON-RPC listener on 127.0.0.1:9090 (wait a few seconds for the window).
It refuses to start if another Coot already holds port 9090.
Tools (from Phenix's bridge): `coot_ping`, `run_python` (single expression; wrap lists/tuples in `str()`),
`run_python_multiline`, `search_coot_functions`, `coot_info`. Notes from Coot's developers are in
`$PHENIX/share/coot/mcp/skills/`.

## Gotchas found while building it (Coot 1.3.4, Phenix 2.3rc2)

- `coot_headless_api` must run under `phenix.python`; plain `bin/python` stalls for minutes in the Coot constructor.
- The Coot GUI needs `COOT_PREFIX` and `COOT_DATA_DIR` set, or monomer dictionaries fail to load (`coot-live` sets them).
- `fit_ligand_right_here` segfaults when a difference map is loaded, so `coot_fit_ligand` with x/y/z uses
  place → random jiggle → rigid-body instead.
- `flip_peptide_using_cid("//A/n")` flips the n-1/n peptide; the tool takes the first residue and offsets.
  Residues with alt confs need `alt_conf`.
- Undo needs `set_make_backups(True)` (done at start), and undo history is per molecule: run trial edits on a
  `coot_copy_model` copy. Mutate + autofit = two undo steps.
- Headless real-space refinement can drag nearby waters away; check them after `coot_refine`.
- R-factors are not exposed by this build's Python bindings and maps are static: re-run phenix.refine for new maps.
- A Coot window started from a terminal may open behind other windows; bring it forward via `run_python`:
  `list(coot_gui_api.application().get_windows())[0].present()`.
- `coot_gui.interesting_things_gui` is GTK3-only in Coot 1.3.4; build button lists with GTK4 (`gi.repository.Gtk`, `Gtk.Window`).

## Licence

MIT; see `LICENSE`.
