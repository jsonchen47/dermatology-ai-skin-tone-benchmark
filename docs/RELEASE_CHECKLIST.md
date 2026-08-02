# Public Release Checklist

- Obtain coauthor approval to release per-image application outputs and derived tables.
- Confirm that the MIT code license and CC BY 4.0 study-output terms are acceptable to the author team and institution.
- Confirm the exact ChatGPT prompt, model/version identifiers, app versions, test dates, devices, and required metadata inputs.
- Confirm that no source images, credentials, personal paths, or protected data are present.
- Run `Rscript analysis/run_all.R` from a fresh clone and review warnings.
- Review all files in `figures/` against the submitted manuscript.
- Create a public GitHub repository and push the verified package.
- Enable the repository in Zenodo, create a GitHub release such as `v1.0.0`, and archive that release.
- Add the version-specific Zenodo DOI to the GitHub release, manuscript Data Sharing Statement, and citation metadata.
- Replace the placeholders in `docs/DATA_SHARING_STATEMENT.md` before submission.
- Cite both source datasets according to their original citation requirements.
