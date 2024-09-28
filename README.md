# Neural cell types and circuits for vocal learning

My post-doctoral research used comparative high-throughput transcriptomics to understand how the brain produces complex, learned behaviors like speech and language. These projects were the first to implement single-cell RNA sequencing in songbirds, establishing a template for the field from molecular to computational components, and the results have broad implications for understanding the genetic toolkits that neurons and circuits use to perform advanced computations.

# `pallium`

Single-cell transcriptomics in the songbird pallium as published in:

Colquitt, B. M.\*, Merullo, D. P.\*, et al. (2021). Cellular transcriptomics reveals evolutionary identities of songbird vocal circuits. *Science*. https://doi.org/10.1126/science.abd9704 \*co-first author

This is a reduced repository to present the R code for the figures from the publication. [See the original repository associated with the publication for more complete information](https://github.com/bradleycolquitt/songbird_cells).

- `01-code-published.R`: R code to produce Seurat objects for downstream analysis.

- `02-code-raw.R`: Complete R code to create the figures seen in Supplemental Figure 10.

# `striatum`

Single-nucleus RNA sequencing in the songbird striatum as published in:

Xiao, L., Merullo, D. P., Cao, M., Co, M., Kulkarni, A., Konopka, G., & Roberts, T. F. (2020). Expression of FoxP2 in the basal ganglia regulates vocal motor sequences in the adult songbird. *Nature Communications*. https://doi.org/10.1038/s41467-021-22918-2

This is a reduced repository to present the R code for the figures from the publication. [See the original repository associated with the publication for more complete information](https://github.com/konopkalab/songbird_areax).

- `01-dataclustering.md`: Markdown file of R code to create Seurat objects from raw data.
  
- `02-code.R`: R code to produce all figures (main and supplementary) in the accompanying study. 
