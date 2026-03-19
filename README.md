# D-PRISM-T2D-PRS
Reference: Huerta-Chagoya A., Kim J., et al., The Lancet of Diabetes and Endocrinology, 2025. [medRxiv preprint](https://www.medrxiv.org/content/10.1101/2025.07.21.25331778v1)

## Overview

**[Diabetes Polygenic Risk Scores in Multiple ancestries (D-PRISM)](https://primedconsortium.org/study/d-prism)** is an international consortium to develop PRS for different types of diabetes and progression across the lifespan in diverse ancestries. 

As part of D-PRISM, this paper addresses a critical gap by delivering the most comprehensive and rigorously tested set of multi-ancestry PRSs for T2D. By making the PRS weights publicly available through the [PGS Catalog](https://www.pgscatalog.org/publication/PGP000773/), we hope to provide a valuable resource for researchers and clinicians seeking to advance genetic risk stratification and develop prevention strategies for T2D.

## Workflow

### Tools

- Plink (>=1.9)
- R (>=3.5)

### 1. Get necessary weights for your target cohort.

These are **multi-ancestry PRSs** trained for five ancestries separately, yet each leverage genetic information from all ancestries. Therefore, before applying them to downstream analyses, you have to combine the **five ancestry-specific scores into a single metascore**.

First, determine which multi-ancestry PRS you want to apply and download the corresponding set of five ancestry-specific weight files required to construct the single metascore for that specific ancestry group.

For example:
  if the closest genetic ancestry of your target cohort is the Admixed American (AMR), you might want to apply the multi-ancestry PRS trained for AMR. Thus, download `PGS005353-5357` files which include the *ancestry-specific weights for AMR-trained PRS*.

See below for a full list:

  
  | Ancestry where the multi-ancestry PRS was trained | PGS IDs | PGS Names |
  | -------------- | -------------- | -------------- | ------------- |
  | African or American (AFR) | PGS005353-5357 | `DPRISM_T2DPRS_trainedforAFR_[AFR/AMR/EAS/EUR/SAS]weights` |
  | Admixed American (AMR) | PGS005358-5362 | `DPRISM_T2DPRS_trainedforAMR_[AFR/AMR/EAS/EUR/SAS]weights` |
  | East Asian (EAS) | PGS005363-5367 | `DPRISM_T2DPRS_trainedforEAS_[AFR/AMR/EAS/EUR/SAS]weights` |
  | European (EUR) | PGS005368-5372 | `DPRISM_T2DPRS_trainedforEUR_[AFR/AMR/EAS/EUR/SAS]weights` |
  | South Asian (SAS) | PGS005373-5377 | `DPRISM_T2DPRS_trainedforSAS_[AFR/AMR/EAS/EUR/SAS]weights` |
  
### 2. Generate individual ancestry-specific scores.

Each PGS file contains metadata lines starting with `#`. You must format the files to remove them and retain columns for variant ID, effect allele and weight. See above an example (adjust based on the the variant ID in your genotype files):

  ```bash
  #Example for AMR-trained score files assuming bim files have variant IDs in the chr:position formatted
  #use position from $8 if your genotypes are in GRCh38 or from $2 if they are in GRCh37
  zcat PGS005358_hmPOS_GRCh38.txt.gz | grep -v "^#" | awk 'NR>=1 {print $1":"$8, $3, $8 }' | gzip > DPRISM_T2DPRS_trainedforAMR_AFRweights.txt.gz
  zcat PGS005359_hmPOS_GRCh38.txt.gz | grep -v "^#" | awk 'NR>=1 {print $1":"$8, $3, $8 }' | gzip > DPRISM_T2DPRS_trainedforAMR_AMRweights.txt.gz
  zcat PGS005360_hmPOS_GRCh38.txt.gz | grep -v "^#" | awk 'NR>=1 {print $1":"$8, $3, $8 }' | gzip > DPRISM_T2DPRS_trainedforAMR_EASweights.txt.gz
  zcat PGS005361_hmPOS_GRCh38.txt.gz | grep -v "^#" | awk 'NR>=1 {print $1":"$8, $3, $8 }' | gzip > DPRISM_T2DPRS_trainedforAMR_EURweights.txt.gz
  zcat PGS005362_hmPOS_GRCh38.txt.gz | grep -v "^#" | awk 'NR>=1 {print $1":"$8, $3, $8 }' | gzip > DPRISM_T2DPRS_trainedforAMR_SASweights.txt.gz
```

Apply the five ancestry-specific formatted files to your genotype data using PLINK. This step will produce five separate polygenic score values per individual.

Example of usage:
  
  ```bash
yourcohortprefix="mycohort"
training_ancestry="AMR" 

WEIGHTS=("AFR" "AMR" "EAS" "EUR" "SAS")
for weight_ancestry in "${WEIGHTS[@]}"; do

WEIGHT="DPRISM_T2DPRS_trainedfor${training_ancestry}_${weight_ancestry}weights.txt.gz"
OUTPUT="${yourcohortprefix}_PRStrainedfor${training_ancestry}_${weight_ancestry}weights"

plink \
--bfile "$BFILE_PATH" \
--score "$WEIGHT" 1 2 3 header \ #variant ID in col 1, effect allele in col 2, and effect weight in col 3
--out "$OUTPUT"

done
```

`--bfile` *(required)*: Path to your genotypes. Basic QC is desirable (*e.g.*, exclude variants with poor imputation quality, palindromic, MAF<0.005).

`--score` *(required)*: Formatted ancestry-specific weight file. Each has one line per scored variant. In the example, the variant ID is read from column 1, the effect allele is read from column 2, and the weight associated with the effect allele is read from the column 3. You can change the variant ID as needed to fit your bim file, as well as the column numbers.

`--out`: Path to save the output file. One line per individual.

It is recommended to inspect the `nopred` files. Only a low percentage of predictive variants is expected to be lost. The performance of the multi-ancestry PRS would be compromised otherwise.

In the AMR target cohort example, the following five files are expected:
  
  - `mycohort_PRStrainedforAMR_AFRweights.profile`
- `mycohort_PRStrainedforAMR_AMRweights.profile`
- `mycohort_PRStrainedforAMR_EASweights.profile`
- `mycohort_PRStrainedforAMR_EURweights.profile`
- `mycohort_PRStrainedforAMR_SASweights.profile`

### 3. Combine the individual ancestry-specific scores into a metascore.

The final step involves combining the five separate individual scores into a single metascore. This is achieved by linearly combining the five individual scaled scores generated by PLINK. For every individual, the single metascore is calculated using the following general formula:
  
  $$metascore = (score_{AFR} \times trained coefficient_{AFR}) + (score_{AMR} \times trained coefficient_{AMR}) + (score_{EAS} \times trained coefficient_{EAS}) + (score_{EUR} \times trained coefficient_{EUR}) + (score_{SAS} \times trained coefficient_{SAS})$$
  
  where:
  
  $score_{i}$ is the one of the five individual scores calculated by PLINK for the i-th ancestry-specific score.

  $trained coefficient_{i}$ is the trained ancestry-specific coefficient assigned to the i-th ancestry-specific score.

We  `dprism_t2dprs` R function to perform this step. This function automatically scales the individual ancestry-specific scores and applies the trained ancestry-specific coefficients to generate the single metascore.

```r
devtools::install_github("ahuertach/DPRISMt2d")
library(DPRISMt2d)
```

Example of usage:
  
  ```r
dprism_t2dprs(
  AFR_profile_path = "path/to/yourcohort_AFR.profile", 
  AMR_profile_path = "path/to/yourcohort_AMR.profile",
  EAS_profile_path = "path/to/yourcohort_EAS.profile",
  EUR_profile_path = "path/to/yourcohort_EUR.profile",
  SAS_profile_path = "path/to/yourcohort_SAS.profile",
  targetcohort_ancestry = "AMR",
  output_file_path = "targetcohort_metascore.txt"
)
```

In the AMR target cohort example:
  
  ```r
dprism_t2dprs(
  "mycohort_PRStrainedforAMR_AFRweights.profile",
  "mycohort_PRStrainedforAMR_AFRweights.profile",
  "mycohort_PRStrainedforAMR_AFRweights.profile",
  "mycohort_PRStrainedforAMR_AFRweights.profile",
  "mycohort_PRStrainedforAMR_AFRweights.profile",
  "AMR",
  "mycohort_dprismt2dprs.txt"
)
```

The resulting file `mycohort_dprismt2dprs.txt` will look like this. `METASCORE` is the raw average score and the `ZMETASCORE` is the normalized average score (scaled mean zero).

```text
FID IID METASCORE ZMETASCORE
ID1 ID1 -0.0000000123 1.543
ID2 ID2 0.0000000056 -0.891
```

### 3. Downstream Analysis

Once you have generated your `mycohort_dprismt2dprs.txt` file, check the downstream_analyses/ folder for example script on how to plot the distributions and run association testing. You will need to first merge phenotype and covariate columns to your PRS file.

Example of usage:

```r
prsperformance(
  data = my_data,
  outcome_col = "T2D",
  covariates = c("sex", "age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
  output_prefix = "myprsperformance"
)
```

It will generate standard performance metrics such as OR/SD, AUC with and without the PRS, Nagelkerke r2, as well as PRS density plots in the overall and by status sample, and ROC curves.

## 📬 Contact & Feedback

**Thank you for using our PRSs!** 

If you have any feedback, encounter unexpected errors, please reach out.


