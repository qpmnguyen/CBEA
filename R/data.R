#' Gingival data set from the Human Microbiome Project
#' @usage data(hmp_gingival)
#' @format A list with two elements
#' \describe{
#'     \item{data}{The microbiome relative abundance data with relevant metadata
#'         obtained from the Human Microbiome Project via the \code{HMP16SData} package (snapshot: 11-15-2021).
#'         The data set is hosted the container of type \code{phyloseq}. Using the
#'         \code{mia} package users can convert it to the \code{TreeSummarizedExperiment}
#'         type.}
#'     \item{set}{Sets of microbes based on their metabolism annotation at the Genera level.
#'         Annotations obtained via Calagaro et al.'s repository on Zenodo
#'         (\url{https://doi.org/10.5281/zenodo.3942108})}
#' }
#' @references Data can be downloaded directly from \url{https://hmpdacc.org/hmp/}
#' @references R interface of the data from \url{https://doi.org/doi:10.18129/B9.bioc.HMP16SData}
#' @references Beghini F, Renson A, Zolnik CP, Geistlinger L, Usyk M, Moody TU, et al.
#'     Tobacco Exposure Associated with Oral Microbiota Oxygen Utilization in the
#'     New York City Health and Nutrition Examination Study. Annals of Epidemiology. 2019;34:18–25.e3. doi:10.1016/j.annepidem.2019.03.005
#' @references  Consortium THMP, Huttenhower C, Gevers D, Knight R, Abubucker S, Badger
#'    JH, et al. Structure, Function and Diversity of the Healthy Human Microbiome.
#'    Nature. 2012;486(7402):207–214. doi:10.1038/nature11234.
#' @references Calgaro M, Romualdi C, Waldron L, Risso D, Vitulo N. Assessment of Statistical
#'    Methods from Single Cell, Bulk RNA-Seq, and Metagenomics Applied to
#'    Microbiome Data. Genome Biology. 2020;21(1):191.
#'    doi:10.1186/s13059-020-02104-1
#' @references Schiffer L, Azhar R, Shepherd L, Ramos M, Geistlinger L, Huttenhower C, et al.
#'    HMP16SData: Efficient Access to the Human Microbiome Project through
#'    Bioconductor. American Journal of Epidemiology. 2019;doi:10.1093/aje/kwz006.
"hmp_gingival"
