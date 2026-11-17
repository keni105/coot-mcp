"""coot-edit: an MCP server for model building and validation with headless Coot.

Runs on Phenix's Python (which ships coot_headless_api 1.3.4 and fastmcp):

    $PHENIX/bin/phenix.python coot_edit_server.py      (e.g. PHENIX=/Applications/phenix-2.3rc2-6237)

phenix.python (the libtbx dispatcher) is required: under the bare bin/python
the Coot constructor misses its environment and stalls for minutes.

Molecules are addressed by Coot's integer molecule index (imol), as returned
by the load tools. Residue/atom selections use Coot CIDs, e.g. "//A/15",
"//A/10-20", "//A/15/OG1". Every editing tool reports which residues the CID
actually matched, so a selection that hit the wrong thing is visible.

Nothing is written to disk except by coot_write_model, which never overwrites
a file unless overwrite=True.
"""
import io
import json
import math
import os
import sys
import threading
import traceback
import contextlib

from fastmcp import FastMCP


mcp = FastMCP("coot_edit")

_lock = threading.RLock()     # coot_headless_api is not thread-safe; FastMCP runs tools in threads
_mc = None
_registry = {}                # imol -> {"kind": "model"|"map"|"ligand", "path": str, "label": str, "diff": bool}
_refinement_map = None


# --------------------------------------------------------------------------
# plumbing
# --------------------------------------------------------------------------

def _detach_stdout():
  """Keep Coot's C++ chatter (written straight to fd 1) off the JSON-RPC pipe."""
  sys.stdout.flush()
  real_fd = os.dup(1)
  devnull = os.open(os.devnull, os.O_WRONLY)
  os.dup2(devnull, 1)
  os.close(devnull)
  sys.stdout = io.TextIOWrapper(os.fdopen(real_fd, "wb", buffering=0),
                                encoding="utf-8", line_buffering=True,
                                write_through=True)


def mc():
  global _mc
  if _mc is None:
    import coot_headless_api
    _mc = coot_headless_api.molecules_container_t(False)
    _mc.set_make_backups(True)      # undo/redo need this
    try:
      _mc.set_use_gemmi(True)
    except Exception:
      pass
  return _mc


class CootError(Exception):
  pass


def _model(imol):
  if not mc().is_valid_model_molecule(imol):
    raise CootError("imol %s is not a loaded model (see coot_inventory)" % imol)
  return imol


def _map(imol):
  if not mc().is_valid_map_molecule(imol):
    raise CootError("imol %s is not a loaded map (see coot_inventory)" % imol)
  return imol


def _refmap(imol_map):
  if imol_map is None:
    if _refinement_map is None:
      raise CootError("no map given and no refinement map set; load maps first")
    return _refinement_map
  return _map(imol_map)


def _selection(imol, cid):
  """Residues a CID matches: [(chain, resno, ins, name, atoms[])]."""
  text = mc().get_molecule_selection_as_json(imol, cid)
  if not text:
    return []
  out = []
  for model in json.loads(text).get("models", []):
    for chain in model.get("chains", []):
      for r in chain.get("residues", []):
        out.append((chain.get("chain_id", ""), r.get("seqnum"), r.get("ins_code", ""),
                    r.get("name", ""), r.get("atoms", [])))
  return out


def _require_match(imol, cid, max_residues=None):
  res = _selection(imol, cid)
  if not res:
    raise CootError("CID %r matched nothing in imol %d" % (cid, imol))
  if max_residues is not None and len(res) > max_residues:
    raise CootError("CID %r matched %d residues (limit %d for this tool): %s ..."
                    % (cid, len(res), max_residues, _labels(res[:5])))
  return res


def _labels(res):
  return ["%s/%s%s %s" % (c, n, i or "", name) for c, n, i, name, _ in res]


def _split_cid(cid):
  """'//A/15' or '//A/15A' -> ('A', 15, ''). One residue only."""
  parts = cid.strip("/").split("/")
  if len(parts) < 2:
    raise CootError("expected a residue CID like //A/15, got %r" % cid)
  chain, rn = parts[0], parts[1]
  num = "".join(ch for ch in rn if ch in "-0123456789")
  ins = rn[len(num):].strip(".")
  try:
    return chain, int(num), ins
  except ValueError:
    raise CootError("could not read a residue number from %r" % cid)


def _rspec(s):
  return "%s/%d%s" % (s.chain_id, s.res_no, s.ins_code or "")


def _aspec(s):
  alt = (":" + s.alt_conf) if s.alt_conf else ""
  return "%s %s%s" % (_rspec(s), s.atom_name.strip(), alt)


def _xyz(imol, cid):
  """{(residue label, atom name, altloc): (x, y, z)} for change detection."""
  return {(lab, a["name"], a.get("altLoc", "")): (a["x"], a["y"], a["z"])
          for lab, r in zip(_labels(_selection(imol, cid)), _selection(imol, cid)) for a in r[4]}


def _rmsd(imol_map):
  r = mc().get_map_rmsd_approx(imol_map)
  return r if r and r > 0 else None


def _residue_fit(imol, res, imol_map=None, imol_diff=None):
  """Per-residue density summary in sigma units at the atom centres."""
  m = mc()
  out = []
  imol_map = imol_map if imol_map is not None else _refinement_map
  imol_diff = imol_diff if imol_diff is not None else _default_diff()
  r2 = _rmsd(imol_map) if imol_map is not None else None
  rd = _rmsd(imol_diff) if imol_diff is not None else None
  for c, n, i, name, atoms in res:
    row = {"residue": "%s/%s%s %s" % (c, n, i or "", name), "n_atoms": len(atoms)}
    if atoms and r2:
      vals = [m.get_density_at_position(imol_map, a["x"], a["y"], a["z"]) / r2 for a in atoms]
      row["2fofc_mean_sigma"] = round(sum(vals) / len(vals), 2)
      worst = min(range(len(vals)), key=lambda k: vals[k])
      row["2fofc_min"] = "%s %.2f" % (atoms[worst]["name"], vals[worst])
    if atoms and rd:
      vals = [m.get_density_at_position(imol_diff, a["x"], a["y"], a["z"]) / rd for a in atoms]
      lo = min(range(len(vals)), key=lambda k: vals[k])
      hi = max(range(len(vals)), key=lambda k: vals[k])
      row["fofc_min"] = "%s %.2f" % (atoms[lo]["name"], vals[lo])
      row["fofc_max"] = "%s %.2f" % (atoms[hi]["name"], vals[hi])
    out.append(row)
  return out


def _default_diff():
  for k, v in _registry.items():
    if v["kind"] == "map" and v.get("diff") and mc().is_valid_map_molecule(k):
      return k
  return None


def _edit_report(imol, cid, status, extra=None, before=None):
  """Standard return for an editing tool: status + what it touched + new fit."""
  res = _selection(imol, cid)
  rep = {"status": status, "imol": imol, "cid": cid, "matched": _labels(res)[:40]}
  if len(res) > 40:
    rep["matched_truncated"] = len(res)
  if res and len(res) <= 30 and _refinement_map is not None:
    rep["fit_after"] = _residue_fit(imol, res)
    if before is not None:
      rep["fit_before"] = before
  if extra:
    rep.update(extra)
  return rep


def _fit_before(imol, cid):
  if _refinement_map is None:
    return None
  res = _selection(imol, cid)
  return _residue_fit(imol, res) if res and len(res) <= 30 else None


def _run(fn):
  """Serialise Coot access and turn exceptions into readable errors."""
  with _lock:
    try:
      return fn()
    except CootError as e:
      return {"error": str(e)}
    except Exception as e:
      return {"error": "%s: %s" % (type(e).__name__, e),
              "trace": traceback.format_exc(limit=3)}


def _vinfo_rows(v, worst_n, ascending=True):
  rows = []
  for ch in v.cviv:
    for r in ch.rviv:
      s = r.residue_spec
      rows.append({"residue": "%s/%d%s" % (s.chain_id, s.res_no, s.ins_code or ""),
                   "value": round(float(r.function_value), 4),
                   "label": r.label})
  rows = [r for r in rows if not math.isnan(r["value"])]
  rows.sort(key=lambda r: r["value"], reverse=not ascending)
  return {"n_residues": len(rows), "worst": rows[:worst_n]}


def _place(p):
  d = {"x": round(p.x, 3), "y": round(p.y, 3), "z": round(p.z, 3)}
  for attr in ("feature_value", "badness", "button_label", "feature_type"):
    if hasattr(p, attr):
      val = getattr(p, attr)
      d[attr] = round(val, 3) if isinstance(val, float) else val
  return d


# --------------------------------------------------------------------------
# loading / inventory
# --------------------------------------------------------------------------

@mcp.tool()
def coot_load_model(path: str) -> dict:
  """Read a PDB or mmCIF model. Returns its imol. The first load also starts Coot
  (about a second)."""
  def f():
    p = os.path.abspath(os.path.expanduser(path))
    if not os.path.isfile(p):
      raise CootError("no such file: %s" % p)
    imol = mc().read_coordinates(p)
    if imol < 0:
      raise CootError("Coot could not read %s" % p)
    _registry[imol] = {"kind": "model", "path": p, "label": os.path.basename(p)}
    m = mc()
    return {"imol": imol, "path": p, "n_atoms": m.get_number_of_atoms(imol),
            "chains": list(m.get_chains_in_model(imol)),
            "residues_without_dictionary": list(m.get_residue_names_with_no_dictionary(imol))}
  return _run(f)


@mcp.tool()
def coot_load_maps(mtz_path: str, f: str = "2FOFCWT", phi: str = "PH2FOFCWT",
                   f_diff: str = "FOFCWT", phi_diff: str = "PHFOFCWT") -> dict:
  """Read 2mFo-DFc and mFo-DFc maps from an MTZ (phenix.refine labels by default)
  and make the 2mFo-DFc map the refinement map. Maps are static: after edits,
  re-run phenix.refine for new maps and R-factors. Pass f_diff="" to skip the
  difference map."""
  def f_():
    global _refinement_map
    p = os.path.abspath(os.path.expanduser(mtz_path))
    if not os.path.isfile(p):
      raise CootError("no such file: %s" % p)
    m = mc()
    i2 = m.read_mtz(p, f, phi, "", False, False)
    if i2 < 0:
      raise CootError("could not read %s/%s from %s" % (f, phi, p))
    _registry[i2] = {"kind": "map", "path": p, "label": "%s %s" % (os.path.basename(p), f), "diff": False}
    out = {"two_fofc": i2, "rmsd_2fofc": _rmsd(i2)}
    idiff = -1
    if f_diff:
      idiff = m.read_mtz(p, f_diff, phi_diff, "", False, True)
      if idiff < 0:
        raise CootError("could not read %s/%s from %s" % (f_diff, phi_diff, p))
      _registry[idiff] = {"kind": "map", "path": p, "label": "%s %s" % (os.path.basename(p), f_diff), "diff": True}
      out["difference"] = idiff
      out["rmsd_fofc"] = _rmsd(idiff)
    m.set_imol_refinement_map(i2)
    _refinement_map = i2
    out["refinement_map"] = i2
    return out
  return _run(f_)


@mcp.tool()
def coot_load_dictionary(cif_path: str, imol: int = -999999) -> dict:
  """Read a ligand restraint CIF. imol=-999999 (default) makes it available to
  every molecule; pass a model imol to scope it to that model."""
  def f():
    p = os.path.abspath(os.path.expanduser(cif_path))
    r = mc().import_cif_dictionary(p, imol)
    return {"status": r, "path": p}
  return _run(f)


@mcp.tool()
def coot_get_monomer(comp_id: str, x: float = None, y: float = None, z: float = None,
                     imol_enc: int = -999999) -> dict:
  """Make a new molecule containing one ligand/monomer (from CCP4 monomer
  library or a loaded dictionary). If x/y/z are given it is placed there."""
  def f():
    m = mc()
    if x is not None:
      imol = m.get_monomer_and_position_at(comp_id, imol_enc, x, y, z)
    else:
      imol = m.get_monomer_from_dictionary(comp_id, imol_enc, True)
    if imol < 0:
      raise CootError("no dictionary for %s" % comp_id)
    _registry[imol] = {"kind": "ligand", "path": "", "label": comp_id}
    return {"imol": imol, "n_atoms": m.get_number_of_atoms(imol)}
  return _run(f)


@mcp.tool()
def coot_inventory() -> dict:
  """Everything loaded, with current atom counts, the refinement map and map weight."""
  def f():
    m = mc() if _registry else None
    items = []
    for k, v in sorted(_registry.items()):
      d = dict(v, imol=k)
      if m and m.is_valid_model_molecule(k):
        d["n_atoms"] = m.get_number_of_atoms(k)
      elif m and m.is_valid_map_molecule(k):
        d["rmsd"] = _rmsd(k)
      else:
        d["closed"] = True
      items.append(d)
    return {"molecules": items, "refinement_map": _refinement_map,
            "map_weight": m.get_map_weight() if m else None}
  return _run(f)


@mcp.tool()
def coot_close(imol: int) -> dict:
  """Close one molecule (model or map)."""
  def f():
    global _refinement_map
    r = mc().close_molecule(imol)
    _registry.pop(imol, None)
    if imol == _refinement_map:
      _refinement_map = None
    return {"status": r}
  return _run(f)


# --------------------------------------------------------------------------
# inspection / validation
# --------------------------------------------------------------------------

@mcp.tool()
def coot_residue_info(imol: int, cid: str, include_atoms: bool = False) -> dict:
  """Residues matched by a CID with per-residue density fit in sigma:
  2mFo-DFc mean and weakest atom, mFo-DFc most negative and most positive atom.
  include_atoms=True also returns atom names, coordinates, B and occupancy."""
  def f():
    _model(imol)
    res = _require_match(imol, cid, max_residues=200)
    out = {"residues": _residue_fit(imol, res)}
    if include_atoms:
      out["atoms"] = {lab: [{k: a.get(k) for k in ("name", "altLoc", "x", "y", "z", "tempFactor", "occupancy", "element") if k in a}
                            for a in r[4]] for lab, r in zip(_labels(res), res)}
    return out
  return _run(f)


@mcp.tool()
def coot_density_at(points: list[list[float]], imol_map: int = None) -> dict:
  """Density in sigma (and raw) at a list of [x, y, z] points. Defaults to the
  refinement map; pass the difference map imol for mFo-DFc."""
  def f():
    im = _refmap(imol_map)
    r = _rmsd(im)
    m = mc()
    return {"imol_map": im, "rmsd": r, "readings": [
      {"xyz": p, "raw": round(m.get_density_at_position(im, *p), 4),
       "sigma": round(m.get_density_at_position(im, *p) / r, 2) if r else None} for p in points]}
  return _run(f)


@mcp.tool()
def coot_validate(imol: int, imol_map: int = None, analyses: list[str] = None,
                  worst_n: int = 15) -> dict:
  """Per-residue validation, worst residues first. analyses is any of:
  "density_fit" (higher = worse), "density_correlation" (lower = worse),
  "rotamer", "ramachandran", "omega", "q_score" (lower = worse).
  Defaults to density_correlation, rotamer, ramachandran."""
  def f():
    _model(imol)
    m = mc()
    im = None
    want = analyses or ["density_correlation", "rotamer", "ramachandran"]
    out = {}
    for a in want:
      if a in ("density_fit", "density_correlation", "q_score"):
        im = im if im is not None else _refmap(imol_map)
      if a == "density_fit":
        out[a] = _vinfo_rows(m.density_fit_analysis(imol, im), worst_n, ascending=False)
      elif a == "density_correlation":
        out[a] = _vinfo_rows(m.density_correlation_analysis(imol, im), worst_n)
      elif a == "q_score":
        out[a] = _vinfo_rows(m.get_q_score(imol, im), worst_n)
      elif a == "rotamer":
        out[a] = _vinfo_rows(m.rotamer_analysis(imol), worst_n)
      elif a == "ramachandran":
        out[a] = _vinfo_rows(m.ramachandran_analysis(imol), worst_n)
      elif a == "omega":
        out[a] = _vinfo_rows(m.peptide_omega_analysis(imol), worst_n, ascending=False)
      else:
        out[a] = {"error": "unknown analysis"}
    return out
  return _run(f)


@mcp.tool()
def coot_geometry_outliers(imol: int, cid: str, include_non_bonded: bool = False,
                           include_hydrogens: bool = False, worst_n: int = 20) -> dict:
  """Restraint-geometry distortions (bonds, angles, planes, chirals, optionally
  clashes) for a selection, against the dictionaries, worst first. Restraints
  involving riding hydrogens are skipped unless include_hydrogens=True."""
  def f():
    _model(imol)
    _require_match(imol, cid)
    vv = mc().get_validation_vs_dictionary_for_selection(imol, cid, include_non_bonded)
    rows = []
    for cont in vv:
      for k in range(cont.size()):
        g = cont.get_geometry_distortion_info(k)
        atoms = [_aspec(a) for a in g.atom_specs]
        if not include_hydrogens and any(a.atom_name.strip().startswith(("H", "D")) and len(a.atom_name.strip()) > 0
                                         and a.atom_name.strip()[0] in "HD" for a in g.atom_specs):
          continue
        rows.append({"distortion": round(g.get_distortion(), 2),
                     "type": str(g.restraint.restraint_type).split(".")[-1],
                     "target": round(g.restraint.target_value, 3),
                     "atoms": atoms})
    rows.sort(key=lambda r: -r.get("distortion", 0))
    return {"n": len(rows), "worst": rows[:worst_n]}
  return _run(f)


@mcp.tool()
def coot_difference_peaks(imol_model: int, imol_diff_map: int = None,
                          n_rmsd: float = 4.0, max_peaks: int = 30) -> dict:
  """Peaks in the mFo-DFc map above n_rmsd (positive and negative), strongest
  first, with height in sigma and the residue Coot associates with each."""
  def f():
    _model(imol_model)
    idm = imol_diff_map if imol_diff_map is not None else _default_diff()
    if idm is None:
      raise CootError("no difference map loaded")
    _map(idm)
    peaks = mc().difference_map_peaks(idm, imol_model, n_rmsd)
    rd = _rmsd(idm)
    rows = []
    for p in peaks:
      d = _place(p)
      d.pop("button_label", None)
      if rd:
        d["sigma"] = round(p.feature_value / rd, 2)
      try:
        d["near"] = _rspec(p.residue_spec)
      except Exception:
        pass
      rows.append(d)
    rows.sort(key=lambda r: -abs(r.get("sigma", r.get("feature_value", 0))))
    return {"n": len(rows), "peaks": rows[:max_peaks]}
  return _run(f)


@mcp.tool()
def coot_unmodelled_blobs(imol_model: int, imol_map: int = None,
                          rmsd_cut_off: float = 1.4, max_blobs: int = 20) -> dict:
  """Clusters of 2mFo-DFc density not covered by the model (candidate ligands,
  ions, missing side chains), biggest first."""
  def f():
    _model(imol_model)
    blobs = mc().unmodelled_blobs(imol_model, _refmap(imol_map), rmsd_cut_off)
    rows = [_place(b) for b in blobs]
    rows.sort(key=lambda r: -r.get("feature_value", 0))
    return {"n": len(rows), "blobs": rows[:max_blobs]}
  return _run(f)


@mcp.tool()
def coot_check_waters(imol_model: int, imol_map: int = None, b_factor_lim: float = 60.0,
                      outlier_sigma_level: float = 0.8, min_dist: float = 2.3,
                      max_dist: float = 3.5) -> dict:
  """Suspicious waters: high B, weak density, or too close/too far from
  H-bond partners (any of these)."""
  def f():
    _model(imol_model)
    bad = mc().find_water_baddies(imol_model, _refmap(imol_map), b_factor_lim,
                                  outlier_sigma_level, min_dist, max_dist, False, False)
    return {"n": len(bad), "waters": [{"water": _aspec(a), "why": getattr(a, "string_user_data", "")} for a in bad]}
  return _run(f)


@mcp.tool()
def coot_clashes(imol: int, ligand_cid: str = "", max_n: int = 30) -> dict:
  """Atom overlaps (clashes) for the whole model, or around one ligand."""
  def f():
    _model(imol)
    m = mc()
    ov = m.get_overlaps_for_ligand(imol, ligand_cid) if ligand_cid else m.get_atom_overlaps(imol)
    rows = []
    for o in ov:
      rows.append({"atom_1": _aspec(o.atom_spec_1), "atom_2": _aspec(o.atom_spec_2),
                   "overlap_volume": round(o.overlap_volume, 3),
                   "h_bond": bool(getattr(o, "is_h_bond", False))})
    rows.sort(key=lambda r: -r["overlap_volume"])
    return {"n": len(rows), "worst": rows[:max_n]}
  return _run(f)


@mcp.tool()
def coot_residues_near(imol: int, cid: str, distance: float = 5.0) -> dict:
  """Residues within `distance` Å of a residue (e.g. coordination shell)."""
  def f():
    _model(imol)
    _require_match(imol, cid, max_residues=1)
    near = mc().get_residues_near_residue(imol, cid, distance)
    return {"residues": [_rspec(r) for r in near]}
  return _run(f)


# --------------------------------------------------------------------------
# editing
# --------------------------------------------------------------------------

@mcp.tool()
def coot_refine(imol: int, cid: str, mode: str = "TRIPLE", n_cycles: int = 1000) -> dict:
  """Real-space refine against the refinement map. mode: SINGLE, TRIPLE,
  QUINTUPLE, HEPTUPLE, SPHERE, BIG_SPHERE, CHAIN, ALL (relative to the residue in
  cid), or LITERAL to refine exactly the residues the cid matches (e.g. a range
  //A/10-20). Reports density fit before and after."""
  def f():
    _model(imol)
    if _refinement_map is None:
      raise CootError("load maps first (coot_load_maps)")
    m = mc()
    before = _fit_before(imol, cid)
    if mode.upper() == "LITERAL":
      _require_match(imol, cid, max_residues=300)
      st = m.refine_residues_using_atom_cid(imol, cid, "LITERAL", n_cycles)
    else:
      _require_match(imol, cid, max_residues=1)
      st = m.refine_residues_using_atom_cid(imol, cid, mode.upper(), n_cycles)
    return _edit_report(imol, cid, st, before=before)
  return _run(f)


@mcp.tool()
def coot_set_refinement_options(map_weight: float = None, use_torsion_restraints: bool = None,
                                use_rama_restraints: bool = None, rama_weight: float = None) -> dict:
  """Change real-space refinement settings; returns the current map weight."""
  def f():
    m = mc()
    if map_weight is not None:
      m.set_map_weight(map_weight)
    if use_torsion_restraints is not None:
      m.set_use_torsion_restraints(use_torsion_restraints)
    if use_rama_restraints is not None:
      m.set_use_rama_plot_restraints(use_rama_restraints)
    if rama_weight is not None:
      m.set_rama_plot_restraints_weight(rama_weight)
    return {"map_weight": m.get_map_weight()}
  return _run(f)


@mcp.tool()
def coot_mutate(imol: int, cid: str, new_residue_type: str, autofit: bool = True) -> dict:
  """Mutate one residue to a 3-letter type (e.g. "TYR"), then (by default) fit
  its rotamer into the refinement map. With autofit, reverting takes two
  coot_undo calls (one for the fit, one for the mutation)."""
  def f():
    _model(imol)
    res = _require_match(imol, cid, max_residues=1)
    old = res[0][3]
    st = mc().mutate(imol, cid, new_residue_type.upper())
    extra = {"from": old, "to": new_residue_type.upper()}
    if st and autofit and _refinement_map is not None:
      c, n, i = _split_cid(cid)
      extra["autofit_status"] = mc().auto_fit_rotamer(imol, c, n, i, "", _refinement_map)
    return _edit_report(imol, cid, st, extra)
  return _run(f)


@mcp.tool()
def coot_auto_fit_rotamer(imol: int, cid: str, alt_conf: str = "") -> dict:
  """Pick and fit the best side-chain rotamer for one residue into the refinement map."""
  def f():
    _model(imol)
    _require_match(imol, cid, max_residues=1)
    before = _fit_before(imol, cid)
    c, n, i = _split_cid(cid)
    st = mc().auto_fit_rotamer(imol, c, n, i, alt_conf, _refmap(None))
    return _edit_report(imol, cid, st, before=before)
  return _run(f)


@mcp.tool()
def coot_flip_peptide(imol: int, cid: str, alt_conf: str = "") -> dict:
  """Flip the peptide after residue cid, e.g. "//A/15" flips the 15-16 peptide
  (carbonyl O of 15 swaps side). Follow with coot_refine."""
  def f():
    _model(imol)
    _require_match(imol, cid, max_residues=1)
    c, n, i = _split_cid(cid)
    pair = "//%s/%d-%d" % (c, n, n + 1)
    before = _fit_before(imol, pair)
    _require_match(imol, "//%s/%d" % (c, n + 1), max_residues=1)
    # Coot names a peptide by its second residue: //A/16 flips the O of 15.
    xyz0 = _xyz(imol, pair)
    st = mc().flip_peptide_using_cid(imol, "//%s/%d%s" % (c, n + 1, i), alt_conf)
    moved = [k for k, v in _xyz(imol, pair).items() if xyz0.get(k) != v]
    extra = {"atoms_moved": moved}
    if not moved:
      alts = sorted({k[2] for k in xyz0 if k[2]})
      extra["hint"] = ("nothing moved; these residues have alt confs %s - pass alt_conf" % alts
                       if alts else "nothing moved (chain break or terminus?)")
    return _edit_report(imol, pair, st, extra, before=before)
  return _run(f)


@mcp.tool()
def coot_side_chain_180(imol: int, atom_cid: str) -> dict:
  """Rotate the last chi by 180° (His/Asn/Gln flips), e.g. "//A/42/ND2"."""
  def f():
    _model(imol)
    rcid = "/".join(atom_cid.split("/")[:4])
    _require_match(imol, rcid, max_residues=1)
    before = _fit_before(imol, rcid)
    st = mc().side_chain_180(imol, atom_cid)
    return _edit_report(imol, rcid, st, before=before)
  return _run(f)


@mcp.tool()
def coot_fill_partial_residues(imol: int) -> dict:
  """Rebuild every residue with missing side-chain atoms and fit them to the map."""
  def f():
    _model(imol)
    m = mc()
    missing = [_rspec(r) for r in m.residues_with_missing_atoms(imol)]
    st = m.fill_partial_residues(imol)
    return {"status": st, "rebuilt": missing}
  return _run(f)


@mcp.tool()
def coot_delete(imol: int, cid: str, scope: str = "RESIDUE") -> dict:
  """Delete by CID. scope: ATOM, RESIDUE, WATER, CHAIN, MOLECULE, LITERAL.
  Reports what was matched before deleting."""
  def f():
    _model(imol)
    res = _require_match(imol, cid)
    n0 = mc().get_number_of_atoms(imol)
    st = mc().delete_using_cid(imol, cid, scope.upper())
    return {"status": list(st) if isinstance(st, tuple) else st, "deleted_from": _labels(res)[:50],
            "n_atoms_before": n0, "n_atoms_after": mc().get_number_of_atoms(imol)}
  return _run(f)


@mcp.tool()
def coot_add_waters(imol_model: int, imol_map: int = None) -> dict:
  """Automatically pick waters into the map. Returns how many were added."""
  def f():
    _model(imol_model)
    n = mc().add_waters(imol_model, _refmap(imol_map))
    return {"waters_added": n}
  return _run(f)


@mcp.tool()
def coot_add_terminal_residue(imol: int, cid: str) -> dict:
  """Build one residue onto the N- or C-terminus at cid (fitted to the map)."""
  def f():
    _model(imol)
    _require_match(imol, cid, max_residues=1)
    st = mc().add_terminal_residue_directly_using_cid(imol, cid)
    return {"status": st, "n_atoms": mc().get_number_of_atoms(imol)}
  return _run(f)


@mcp.tool()
def coot_add_alt_conf(imol: int, cid: str, use_difference_map: bool = True) -> dict:
  """Split a residue into A/B conformers (guided by the difference map if loaded)."""
  def f():
    _model(imol)
    _require_match(imol, cid, max_residues=1)
    idm = _default_diff()
    n0 = mc().get_number_of_atoms_in_residue(imol, cid)
    if use_difference_map and idm is not None:
      mc().split_residue_using_map(imol, cid, idm)
    else:
      mc().add_alternative_conformation(imol, cid)
    n1 = mc().get_number_of_atoms_in_residue(imol, cid)
    return _edit_report(imol, cid, int(n1 > n0), {"atoms_in_residue_before": n0, "atoms_in_residue_after": n1,
                        "note": "if before == after the residue already had alt confs or the split failed"})
  return _run(f)


@mcp.tool()
def coot_rigid_body_fit(imol: int, cid: str, imol_map: int = None) -> dict:
  """Rigid-body fit a selection (residues, a domain "//A/1-120", or a ligand)."""
  def f():
    _model(imol)
    _require_match(imol, cid)
    before = _fit_before(imol, cid)
    st = mc().rigid_body_fit(imol, cid, _refmap(imol_map))
    return _edit_report(imol, cid, st, before=before)
  return _run(f)


@mcp.tool()
def coot_set_occupancy_b(imol: int, cid: str, occupancy: float = None, b_factor: float = None) -> dict:
  """Set occupancy and/or B-factor for every atom a CID matches."""
  def f():
    _model(imol)
    res = _require_match(imol, cid)
    if occupancy is not None:
      mc().set_occupancy(imol, cid, occupancy)
    if b_factor is not None:
      mc().set_temperature_factors_using_cid(imol, cid, b_factor)
    return {"changed": _labels(res)[:50]}
  return _run(f)


def _heavy_sigma(imol, cid="//"):
  r = _rmsd(_refinement_map)
  at = [a for res in _selection(imol, cid) for a in res[4] if a.get("element") not in ("H", "D")]
  if not at or not r:
    return None
  return round(sum(mc().get_density_at_position(_refinement_map, a["x"], a["y"], a["z"]) for a in at) / len(at) / r, 2)


@mcp.tool()
def coot_fit_ligand(imol_ligand: int, x: float = None, y: float = None, z: float = None,
                    imol_protein: int = -1, n_trials: int = 300, n_rounds: int = 3,
                    n_rmsd: float = 1.0, use_conformers: bool = False, n_conformers: int = 10) -> dict:
  """Fit a ligand molecule (from coot_get_monomer or a loaded file) into the
  refinement map.

  With x/y/z: moves the ligand to that point and does rounds of random-jiggle
  + rigid-body fitting, keeping imol_ligand as the result (Coot's own
  fit_ligand_right_here segfaults in this build when a difference map is
  loaded, so it is not used). Reports mean 2mFo-DFc sigma over heavy atoms per
  round - compare with ~1.5-3 sigma for a well-placed ligand.

  Without x/y/z: Coot's whole-map ligand search (needs imol_protein); each
  solution is a new molecule. Merge the one you want with coot_merge, then
  coot_refine it."""
  def f():
    _model(imol_ligand)
    m = mc()
    if _refinement_map is None:
      raise CootError("load maps first (coot_load_maps)")
    if x is not None:
      m.move_molecule_to_new_centre(imol_ligand, x, y, z)
      scores = [{"round": 0, "sigma": _heavy_sigma(imol_ligand)}]
      for k in range(max(1, n_rounds)):
        m.fit_to_map_by_random_jiggle_using_cid(imol_ligand, "//", n_trials, 1.0 if k == 0 else 0.5)
        m.rigid_body_fit(imol_ligand, "//", _refinement_map)
        scores.append({"round": k + 1, "sigma": _heavy_sigma(imol_ligand)})
      return {"imol": imol_ligand, "scores": scores, "fit": _residue_fit(imol_ligand, _selection(imol_ligand, "//"))}
    _model(imol_protein)
    sols = m.fit_ligand(imol_protein, _refinement_map, imol_ligand, n_rmsd, use_conformers, n_conformers)
    out = []
    for s_ in sols:
      if m.get_number_of_atoms(s_.imol) == 0:
        continue
      _registry[s_.imol] = {"kind": "ligand", "path": "", "label": "fit of %d" % imol_ligand}
      out.append({"imol": s_.imol, "sigma": _heavy_sigma(s_.imol),
                  "residues": _labels(_selection(s_.imol, "//"))})
    out.sort(key=lambda d: -(d["sigma"] or 0))
    return {"solutions": out}
  return _run(f)


@mcp.tool()
def coot_merge(imol: int, others: list[int]) -> dict:
  """Merge other molecules (ligands, waters, fragments) into model imol."""
  def f():
    _model(imol)
    st, info = mc().merge_molecules(imol, ":".join(str(o) for o in others))
    return {"status": st, "added": [_rspec(i.spec) for i in info],
            "n_atoms": mc().get_number_of_atoms(imol)}
  return _run(f)


@mcp.tool()
def coot_undo(imol: int, redo: bool = False) -> dict:
  """Undo (or redo) the last edit on a model."""
  def f():
    _model(imol)
    n0 = mc().get_number_of_atoms(imol)
    mc().redo(imol) if redo else mc().undo(imol)
    return {"n_atoms_before": n0, "n_atoms_after": mc().get_number_of_atoms(imol),
            "note": "Coot's return code is not meaningful here; check the residue with coot_residue_info"}
  return _run(f)


@mcp.tool()
def coot_copy_model(imol: int) -> dict:
  """Duplicate a model (a cheap checkpoint before risky edits). Returns the new imol."""
  def f():
    _model(imol)
    new = mc().copy_molecule(imol)
    _registry[new] = dict(_registry.get(imol, {"kind": "model", "path": ""}),
                          label="copy of %d" % imol)
    return {"imol": new}
  return _run(f)


@mcp.tool()
def coot_write_model(imol: int, path: str, overwrite: bool = False) -> dict:
  """Write a model to PDB or mmCIF (by extension). Refuses to overwrite an
  existing file unless overwrite=True."""
  def f():
    _model(imol)
    p = os.path.abspath(os.path.expanduser(path))
    if os.path.exists(p) and not overwrite:
      raise CootError("%s exists; pass overwrite=True or choose a new name" % p)
    mc().write_coordinates(imol, p)
    if not os.path.exists(p):
      raise CootError("Coot did not write %s" % p)
    return {"written": p, "bytes": os.path.getsize(p), "n_atoms": mc().get_number_of_atoms(imol)}
  return _run(f)


# --------------------------------------------------------------------------
# escape hatch
# --------------------------------------------------------------------------

@mcp.tool()
def coot_run_python(code: str) -> dict:
  """Run Python against the headless Coot API for anything the other tools
  don't cover. `mc` is the molecules_container_t, `coot_headless_api` is
  imported. Assign to `result` to return a value; print() output is returned.
  Look up functions with coot_api_search first."""
  def f():
    import coot_headless_api
    buf = io.StringIO()
    ns = {"mc": mc(), "coot_headless_api": coot_headless_api, "result": None}
    with contextlib.redirect_stdout(buf):
      exec(code, ns)
    r = ns.get("result")
    try:
      json.dumps(r)
    except TypeError:
      r = repr(r)
    return {"result": r, "stdout": buf.getvalue()[-20000:]}
  return _run(f)


@mcp.tool()
def coot_api_search(pattern: str, max_results: int = 25) -> dict:
  """Search molecules_container_t methods by name/docstring. Space-separated
  words must all match (case-insensitive)."""
  def f():
    import coot_headless_api
    cls = coot_headless_api.molecules_container_t
    words = pattern.lower().split()
    hits = []
    for name in dir(cls):
      if name.startswith("_"):
        continue
      doc = (getattr(cls, name).__doc__ or "")
      if all(w in (name + " " + doc).lower() for w in words):
        hits.append({"name": name, "doc": doc.strip()[:600]})
    return {"n": len(hits), "results": hits[:max_results]}
  return _run(f)


if __name__ == "__main__":
  _detach_stdout()
  mcp.run(show_banner=False)
