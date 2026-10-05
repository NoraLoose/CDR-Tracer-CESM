# Surface pH and dALK diagnostic pipeline

## Dependency order

```
compute_truth_surface.py        # 1. Download truth fields from S3
        │
        ▼
compute_dalk_timeseries.py      # 2a. Compute dALK scalar timeseries (truth + antitracer)
compute_max_dph.py   ◄──┐      # 3.  Compute max(δpH) maps
                         │
compute_ph.py           ─┘      # 2b. Compute CDR-tracer pH 2D fields (independent of 2a)
```

Steps 2a and 2b are independent of each other and can run in parallel.
Step 3 depends on both 2a (for truth `surface.nc`) and 2b (for `ph_{pid}.nc`).

---

## Scripts and outputs

### 1. `compute_truth_surface.py`
Downloads surface fields from the public S3 archive for all OAE truth experiments.

**Reads:** `s3://us-west-2.opendata.source.coop/cworthy/oae-efficiency-atlas/data/`  
**Writes** (per polygon, into truth analysis subdir):
```
surface.nc   —  ALK, ALK_ALT_CO2, PH, PH_ALT_CO2  (time, nlat, nlon)
```

---

### 2a. `compute_dalk_timeseries.py`
Computes surface δALK area-integral timeseries for truth and CDR-tracer.

**Reads:** `surface.nc` (truth, local) + antitracer history files (local)  
**Writes:**
```
{truth_subdir}/dalk_timeseries.nc        —  surf_true  (time,)  [mmol m⁻¹]
{antitracer_subdir}/dalk_timeseries_{pid}.nc  —  surf_anti  (time,)  [mmol m⁻¹]
```

---

### 2b. `compute_ph.py`
Computes CDR-tracer surface pH 2D fields using PyCO2SYS.

**Reads:** control monthly + nday1 files (local) + antitracer history files (local)  
**Writes** (per polygon, into antitracer subdir):
```
ph_{pid}.nc   —  ph_ctrl, ph_cdr  (time, nlat, nlon)
```

---

### 3. `compute_max_dph.py`
Computes per-cell temporal maximum of the pH perturbation.

**Reads:** `surface.nc` (truth) + `ph_{pid}.nc` (CDR tracer)  
**Writes:**
```
{truth_subdir}/dph_max.nc            —  dph_max, dph_max_month  (nlat, nlon)
{antitracer_subdir}/dph_max_{pid}.nc —  dph_max, dph_max_month  (nlat, nlon)
```

---

## Analysis directory layout

All outputs land under:
```
/global/cfs/projectdirs/m4746/Users/nora/Ocean-CDR-Atlas-v0/data/analysis/
```

Per-polygon truth subdirs:
```
smyle.cdr-atlas-v0.glb-oae_{basin}_{p}_1999-01-01_{mult}.001/
    surface.nc
    dalk_timeseries.nc
    dph_max.nc
```

Shared antitracer subdir:
```
smyle.cdr-atlas-v0.glb-antitracer_1999-01_all-oae-daily.001/
    dalk_timeseries_{pid}.nc
    ph_{pid}.nc
    dph_max_{pid}.nc
```
