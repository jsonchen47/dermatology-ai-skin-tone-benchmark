# Public Release Checklist

- Obtain coauthor approval to release per-image application outputs and derived tables.
- Confirm that the MIT code license and CC BY 4.0 study-output terms are acceptable to the author team and institution.
- Confirm the exact ChatGPT prompt, model/version identifiers, app versions, test dates, devices, and required metadata inputs.
- Confirm that no source images, credentials, personal paths, or protected data are present.
- [x] Run `Rscript analysis/run_all.R` and confirm successful completion.
- Review all files in `figures/` against the submitted manuscript.
- [x] Create a public GitHub repository at https://github.com/jsonchen47/dermatology-ai-skin-tone-benchmark.
- [x] Enable the repository in Zenodo.
- Create the GitHub `v1.0.0` release and confirm that Zenodo archives it.
- Add the version-specific Zenodo DOI to the GitHub release, manuscript Data Sharing Statement, and citation metadata.
- Add the version-specific Zenodo DOI to `docs/DATA_SHARING_STATEMENT.md` before submission.
- Cite both source datasets according to their original citation requirements.
