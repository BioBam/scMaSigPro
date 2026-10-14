# Disabled GitHub Actions

Every workflow in this directory is switched off: GitHub Actions only reads
`.github/workflows/`, so nothing here runs on any push, pull request,
schedule, release, comment or manual dispatch. `dependabot.yml` is parked
here for the same reason.

The files are unchanged. To turn CI back on, move the directory back:

```sh
git mv .github/workflows-disabled .github/workflows
git mv .github/workflows/dependabot.yml .github/dependabot.yml
git rm .github/workflows/README.md
```
