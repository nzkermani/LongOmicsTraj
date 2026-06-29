# LongOmicsTraj (LOT)

**Topology-based analysis of longitudinal omics trajectories**

---

## 🧠 Why LongOmicsTraj?

![Why LongOmicsTraj](figures/why_longomicstraj.png)

---

## 🔬 Core idea

Most methods answer:

> *“Does expression change?”* (Mixed models)
> *“What shape does it follow?”* (GAMs)
> *“Which genes are similar?”* (K-means)

**LOT answers:**

> **“Which genes share the same trajectory topology?”**

---
<p align="center">
  <img src="figures/LongOmicsTraj_overview.png" width="1000">
</p>

## ⚙️ Workflow

```r
library(LongOmicsTraj)

data("lot_glucold")

visits <- c("Baseline","6 Months","30 Months")

lot <- lot_from_empirical(lot_glucold, "transcriptomics", visits, overwrite=TRUE)
lot <- lot_compute_deltas(lot, "transcriptomics", visits, overwrite=TRUE)
lot <- lot_compute_thresholds(lot, "transcriptomics", visits=visits, overwrite=TRUE)
lot <- lot_map_to_ots(lot, visits=visits, overwrite=TRUE)
```

---

## 📊 Output: Trajectory Topology (OTS)

Each gene is mapped to a **temporal behaviour pattern**:

* `down_up` → rebound
* `up_down` → transient
* `flat_flat` → stable
* `down_down` → sustained suppression
* `up_up` → sustained activation

---

## 📈 Example result

![OTS distribution](figures/ots_distribution.png)

LOT reveals:

* Stable programs (`flat_flat`)
* Rebound dynamics (`down_up`)
* Transient responses (`up_down`)

---

## 🚀 What LOT adds

| Method       | What it captures       | What it misses   |
| ------------ | ---------------------- | ---------------- |
| Mixed models | statistical change     | trajectory shape |
| GAMs         | smooth patterns        | shared structure |
| K-means      | similarity             | topology         |
| **LOT**      | **topology of change** | —                |

---

## 🧬 Key insight

> LOT transforms time into **interpretable biological behaviour**.

---

## 📦 Installation

```r
devtools::install_github("nzkermani/LongOmicsTraj")
```

---

## 📍 Status

* ✔ Stable pipeline
* ✔ Tested (testthat)
* ✔ Reproducible
* 🚧 Active development

---

## 👩‍🔬 Author

Nazanin Zounemat-Kermani
Imperial College London




...
