# Where did we get the data from?

- Rossi et al. 2021
    - Chip-exo data
    - Data was taken from the paper's Supplementary Data 2, sheet "All_sectors_x_all_targets" where we extracted the relevant columns for Fkh1,2

- Lupo et al. 2021
    - ChEC-Seq data
    - Data was obtained from direct email correspondence with the authors. 
        - 1 file with all scores for all genes and one where the top 100 targets were selected
        - cer = Saccharomyces cerevisiae (strain BY4743) and, par = Saccharomyces paradoxus (strain CBS432)
        - Hyc and Hyp refer to Hybrid S. cerevisiae genome and S. paradoxus genome respectively

- Hackett et al. 2020
    - Timecourses of TF induction (overexpression) data (microarray)
    - https://idea.research.calicolabs.com/data
        - We downloaded "Gene expression data in wide format" (should be the same as datasetEV1)
        - We used datasetEV1 because some numbers in the other file did not seem right (E.g. entry LV11 and LV18)
        - This contains the log2 expression response after overexpression
        - We extracted the relevant columns of Fkh1,2

- Data from Kemmeren et al. 2014
    - Deletion experiments
    - https://deleteome.holstegelab.nl/
        - Fkh1 is in the [deleteome_responsive_mutants_ex_wt_var_controls.txt] dataset but Fkh2 is not
            - Fkh2 is not considered a responsive mutant
        - We extracted the Fkh data from the [deleteome_all_mutants_controls.txt] file
        - Expression levels (A), ratios (M), and p values
            - M ratios are log2(deletion mutant / wild type) 


- SGD / GEMMER
    - Venters data on SGD appears to be only those targets that score >= 1 in the 37/25C ratio.
        - so only heat shock upregulation of binding with at least a two-fold difference
        - This makes no sense for our use here so we focus on the UTmax (max of TSS and UAS1 data) data in both 25C and 37C taking only those above the 5% FDR cutoff. 
    - For the Mondeel et al. 2019 data SGD contains separate entries for exponential and stationary phase