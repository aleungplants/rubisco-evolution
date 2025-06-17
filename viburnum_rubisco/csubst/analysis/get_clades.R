library(ape)
library(ggtree)
library(dplyr)
library(TreeTools)

here::i_am("analysis/get_clades.R")

tree <- read.tree(here::here("radseq_pruned_csubst.tre")) %>%
  Preorder()

node_labels <- tibble(
  Clade = c("succotinus",
    "lobata",
    "sambucina",
    "coriacea",
    "opulus",
    "oreinotinus",
    "dentata",
    "mollotinus",
    "tinus",
    "solenotinus",
    "lutescentia",
    "euviburnum",
    "lentago",
    "punctata",
    "pseudotinus",
    "urceolata",
    "laminotinus", # higher level clades here and below
    "porphyrotinus",
    "crenotinus",
    "valvotinus"
    ),
  Node = c(getMRCA(tree, c("erosum", "foetidum")),
    getMRCA(tree, c("orientale", "kansuense")),
    getMRCA(tree, c("sambucinum", "beccarii")),
    getMRCA(tree, c("hebanthum", "cylindricu")),
    getMRCA(tree, c("sargentii", "edule")),
    getMRCA(tree, c("undulatum", "loeseneri")),
    getMRCA(tree, c("recognitum", "scabrellum")),
    getMRCA(tree, c("australe", "rafinesqui")),
    getMRCA(tree, c("calvum", "treleasei")),
    getMRCA(tree, c("oliganthum", "sieboldii")),
    getMRCA(tree, c("plicatum", "amplifoliu")),
    getMRCA(tree, c("chinshanen", "cotinifoli")),
    getMRCA(tree, c("rufidulum", "cassinoide")),
    getMRCA(tree, c("punctatum", "lepidotulu")),
    getMRCA(tree, c("sympodiale", "nervosum")),
    getMRCA(tree, c("urceolatum", "taiwanianu")),
    getMRCA(tree, c("erosum", "edule")),
    getMRCA(tree, c("undulatum", "australe")),
    getMRCA(tree, c("oliganthum", "amplifoliu")),
    getMRCA(tree, c("chinshanen", "lepidotulu"))
    )
  ) %>%
  rowwise() %>%
  mutate(Taxa = list(Subtree(tree, Node)$tip.label)) %>%
  tidyr::unnest(Taxa)

node_labels %>%
  readr::write_csv(file = here::here("analysis/landis_clades.csv"))
