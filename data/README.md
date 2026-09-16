# Data

Not included in this repository - this pipeline runs on patient-derived
clinical whole-genome sequencing data from the ICONIK paediatric
cardiomyopathy cohort (trio CRAMs, pedigree, and diagnosis records).

Data is available upon request from the study team.

To run this pipeline yourself, you'll need to supply, in this directory:

- `trios_all37.txt` - tab-separated trio manifest (family ID, then per-role
  individual ID / sample ID / CRAM path, for proband/mother/father)
- `iconik_trioPedigree.csv` - pedigree file with at least `Individual ID`,
  `Sex`, and `Affected Status` columns
- `DRAGEN_complete_trios_plus_siblings.csv` - clinical diagnosis lookup
  (`FamilyNo`, `ICONIK_Family_member`, `ICONIK_proband_diagnosis`)
- `families_data/` - per-family CRAM/CRAI symlinks
- `R135/R135_panel.tsv` and `R135/R135_green_genes.tsv` - the green-gene
  panel, downloadable from [PanelApp](https://panelapp.genomicsengland.co.uk/)
  panel 749

See the scripts in `scripts/SV_script/` for the exact expected file formats
(each script's `--help`/usage comment documents its inputs).
