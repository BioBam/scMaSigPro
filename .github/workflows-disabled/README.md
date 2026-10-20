# Disabled GitHub Actions

Every workflow in this directory is switched off: GitHub Actions only reads
`.github/workflows/`, so nothing here runs on any push, pull request,
schedule, release, comment or manual dispatch. `dependabot.yml` is parked
here for the same reason.

The one live workflow is `.github/workflows/pkgdown.yaml`, which builds the
documentation site and deploys it to the `gh-pages` branch on pushes to
`main`, on published releases and on manual dispatch (maintainer decision,
2026-10-20). Everything else stays here until the maintainer says otherwise.

The parked files are unchanged apart from their branch triggers, which now
name `main` (the `devel` branch was removed on 2026-10-20). To turn the rest
of CI back on:

```sh
git mv .github/workflows-disabled/*.yaml .github/workflows/
git mv .github/workflows-disabled/*.yml .github/workflows/
git mv .github/workflows/dependabot.yml .github/dependabot.yml
git rm .github/workflows-disabled/README.md
```
