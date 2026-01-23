# SpaTalk <img src='https://github.com/ZJUFanLab/SpaTalk/blob/main/img/SpaTalk.png' align="right" height="139" />

<!-- badges: start -->
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/SpaTalk)](https://zaoqu-liu.r-universe.dev/SpaTalk)
[![R ≥ 4.0](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue.svg)](https://cran.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41467--022--32111--8-green)](https://doi.org/10.1038/s41467-022-32111-8)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.6809147.svg)](https://doi.org/10.5281/zenodo.6809147)
<!-- badges: end -->

## Knowledge-Graph-Based Cell-Cell Communication Inference for Spatially Resolved Transcriptomic Data

> **Note**: This repository is maintained by [Zaoqu Liu](https://github.com/Zaoqu-Liu). For the original version, please visit [ZJUFanLab/SpaTalk](https://github.com/ZJUFanLab/SpaTalk).

## Overview

**SpaTalk** is a computational framework for inferring spatially resolved cell-cell communications (CCIs) from spatial transcriptomics (ST) data. The method integrates graph network modeling and knowledge graph approaches to reconstruct ligand-receptor-target signaling networks between spatially proximal cells.

### Key Methodological Features

- **Cell-type Deconvolution**: Non-negative linear model (NNLM) for decomposing spot-based ST data into single-cell resolution
- **Spatial Mapping**: Integration of scRNA-seq reference data with spatial coordinates
- **Graph-based CCI Inference**: Knowledge graph modeling of ligand-receptor-downstream pathway interactions
- **Statistical Validation**: Permutation-based significance testing for identified communications

### Supported Data Types

| Platform | Resolution | Examples |
|----------|-----------|----------|
| Single-cell ST | Cellular | STARmap, MERFISH, seqFISH+ |
| Spot-based ST | Multi-cellular | 10x Visium, Slide-seq, ST |

## Installation

### From R-universe (Recommended)

```r
install.packages("SpaTalk", repos = "https://zaoqu-liu.r-universe.dev")
```

### From GitHub

```r
# Install dependencies
install.packages("devtools")
devtools::install_github("linxihui/NNLM")

# Install SpaTalk
devtools::install_github("Zaoqu-Liu/SpaTalk")
```

### System Requirements

- R ≥ 4.0.0
- Platform: Windows, macOS, Linux
- Dependencies: Seurat (≥3.0.0), NNLM, Matrix, Rcpp

## Quick Start

### 1. Create SpaTalk Object

```r
library(SpaTalk)

# For spot-based ST data
obj <- createSpaTalk(
  st_data = st_counts,      # Gene × Spot count matrix
  st_meta = st_coordinates, # Data frame with 'spot', 'x', 'y'
  species = "Human",        # "Human" or "Mouse"
  if_st_is_sc = FALSE,
  spot_max_cell = 10        # Expected cells per spot
)
```

### 2. Cell-type Deconvolution

```r
# Using built-in NNLM deconvolution
obj <- dec_celltype(
  obj, 
  sc_data = sc_counts,       # scRNA-seq reference

  sc_celltype = cell_labels  # Cell type annotations
)
```

### 3. Infer Cell-Cell Communications

```r
# Load curated databases
data(lrpairs)   # Ligand-receptor pairs from CellTalkDB
data(pathways)  # KEGG/Reactome pathways + AnimalTFDB

# Filter LR pairs with downstream targets
obj <- find_lr_path(obj, lrpairs, pathways)

# Infer CCIs between cell types
obj <- dec_cci(obj, 
  celltype_sender = "Macrophage",
  celltype_receiver = "Fibroblast"
)

# Or infer all pairwise CCIs
obj <- dec_cci_all(obj)
```

### 4. Access Results

```r
# Significant LR pairs with scores
head(obj@lrpair)

# Downstream TF activity
head(obj@tf)
```

## Integrated Databases

SpaTalk incorporates curated biological knowledge from:

| Database | Content | Species |
|----------|---------|---------|
| [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/) | 3,398 human / 2,033 mouse LR pairs | Human, Mouse |
| [KEGG](https://www.kegg.jp/kegg/pathway.html) | Signaling pathways | Human, Mouse |
| [Reactome](https://reactome.org/) | Pathway interactions | Human, Mouse |
| [AnimalTFDB](http://bioinfo.life.hust.edu.cn/AnimalTFDB/) | Transcription factors | Human, Mouse |

## Documentation

- **Tutorial**: [Comprehensive vignette](https://raw.githack.com/multitalk/awesome-cell-cell-communication/main/method/tutorial.html)
- **Wiki**: [Detailed documentation](https://github.com/ZJUFanLab/SpaTalk/wiki)
- **API Reference**: [Function documentation](https://raw.githack.com/ZJUFanLab/SpaTalk/main/vignettes/SpaTalk.pdf)

### Tutorials by Data Type

- [Single-cell ST analysis](https://raw.githack.com/multitalk/awesome-cell-cell-communication/main/method/sc_tutorial.html)
- [Spot-based ST analysis](https://raw.githack.com/multitalk/awesome-cell-cell-communication/main/method/spot_tutorial.html)

## Advanced Features

- **Custom databases**: [Use custom LR pairs](https://github.com/ZJUFanLab/SpaTalk/wiki/Use-customed-lrpairs) | [Use custom pathways](https://github.com/ZJUFanLab/SpaTalk/wiki/Use-customed-pathways)
- **Alternative deconvolution**: [RCTD, Seurat, SPOTlight, stereoscope, cell2location](https://github.com/ZJUFanLab/SpaTalk/wiki/Use-other-deconvolution-methods)
- **Direct inference**: [Skip deconvolution for single-cell ST](https://github.com/ZJUFanLab/SpaTalk/wiki/Directly-infer-cell-cell-communication-skiping-deconvolution)

## Citation

If you use SpaTalk in your research, please cite:

> Shao, X., Li, C., Yang, H., Lu, X., Liao, J., Qian, J., Wang, K., Cheng, J., Yang, P., Chen, H., Xu, X., & Fan, X. (2022). **Knowledge-graph-based cell-cell communication inference for spatially resolved transcriptomic data with SpaTalk.** *Nature Communications*, 13, 4429. https://doi.org/10.1038/s41467-022-32111-8

```bibtex
@article{shao2022spatalk,
  title={Knowledge-graph-based cell-cell communication inference for spatially resolved transcriptomic data with SpaTalk},
  author={Shao, Xin and Li, Chengyu and Yang, Haihong and Lu, Xiaoyan and Liao, Jie and Qian, Jingyang and Wang, Kai and Cheng, Junyun and Yang, Penghui and Chen, Huajun and Xu, Xiao and Fan, Xiaohui},
  journal={Nature Communications},
  volume={13},
  pages={4429},
  year={2022},
  doi={10.1038/s41467-022-32111-8}
}
```

## Maintainers

- **Current maintainer**: [Zaoqu Liu](mailto:liuzaoqu@163.com) (GitHub: [@Zaoqu-Liu](https://github.com/Zaoqu-Liu))
- **Original developer**: [Xin Shao](mailto:xin_shao@zju.edu.cn) (GitHub: [@ZJUFanLab](https://github.com/ZJUFanLab))

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.

---

<p align="center">
  <a href="https://github.com/Zaoqu-Liu/SpaTalk">
    <img src="https://img.shields.io/github/stars/Zaoqu-Liu/SpaTalk?style=social" alt="GitHub stars">
  </a>
</p>
