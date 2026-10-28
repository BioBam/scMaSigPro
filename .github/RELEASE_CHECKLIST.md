# Release checklist

Releases are made from `main` of spriyansh/scMaSigPro and then synced to
BioBam/scMaSigPro. Replace `x.y.z` with the new version.

## Prepare (one pull request to `main`)

- [ ] `Version: x.y.z` in `DESCRIPTION`.
- [ ] `NEWS.md`: the heading reads `# scMaSigPro x.y.z (YYYY-MM-DD)` and
      every change since the last release has its bullet.
- [ ] `data/scmp.ob.RData` regenerated, since the object records the package
      version:
      `Rscript -e 'pkgload::load_all(); source("inst/scripts/make-scmp-ob.R")'`.
- [ ] `inst/CITATION` and the installation line of `README.Rmd` are current;
      `README.md` knitted from `README.Rmd`.
- [ ] `devtools::document()` and `styler::style_pkg(indent_by = 4)` change
      nothing.
- [ ] `devtools::test()` passes and `spelling::spell_check_package()` finds
      nothing.
- [ ] `R CMD build .` and `R CMD check --as-cran scMaSigPro_x.y.z.tar.gz`
      give no ERROR or WARNING and no new NOTE.
- [ ] `BiocCheck::BiocCheck("scMaSigPro_x.y.z.tar.gz")` gives no new ERROR,
      WARNING or NOTE.
- [ ] Every check of the pull request passes, the website build included.

## Release

- [ ] Squash-merge the pull request into `main`.
- [ ] Tag the merge commit and push the tag:
      `git tag -a vx.y.z -m "scMaSigPro x.y.z"` and
      `git push origin vx.y.z`.
- [ ] Publish a GitHub release from the tag on spriyansh/scMaSigPro, with the
      NEWS section of the version as its notes; the website is rebuilt.

## Sync to BioBam/scMaSigPro

- [ ] Open a pull request from `spriyansh:main` to `BioBam:main` and merge it
      with a merge commit, never a squash, so that both repositories keep the
      same commits. BioBam requires signed commits.
- [ ] Push the tag, and any new branch to keep, to BioBam:
      `git push upstream vx.y.z`.
- [ ] Publish the same GitHub release on BioBam/scMaSigPro.
- [ ] On spriyansh/scMaSigPro, use "Sync fork" so that `main` also holds the
      merge commit.

## Later

- [ ] R-universe: once the package is registered, check its build page after
      the release.
- [ ] Bioconductor: from the submission on, versions are `0.99.z` and follow
      the Bioconductor release schedule.
