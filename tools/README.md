# Package verification and Git preparation

[Project home](../README.md) · [Documentation](../docs/README.md) · [Publication instructions](../PUBLICATION.md)

These tools prepare and check the publication package. Scientific reproduction
commands are documented separately in [REPRODUCING.md](../REPRODUCING.md).

| Tool | What it does |
| --- | --- |
| [verify_package.py](verify_package.py) | Read-only verification of the manifest, file inventory, SHA-256 hashes, paths, file sizes, and known credential patterns; optional deeper privacy checks |
| [initialize_anonymous_git.py](initialize_anonymous_git.py) | Verifies an unchanged package, then creates one local commit with a neutral author, email, and date; refuses an existing repository |
| [build_navigation.py](build_navigation.py) | Regenerates the paginated file catalog and readable reference index without network access |

Run from the package root before executing research scripts that regenerate files:

```sh
python3 tools/verify_package.py
```

For deeper privacy verification, provide the private identifiers to check in the
local `ANON_PRIVATE_TOKENS` environment variable as a JSON list of strings, each
at least three characters long. Include relevant names, usernames, hostnames,
email addresses, and other known identifying strings. Keep those values out of
the repository and published reports. Then run:

```sh
python3 tools/verify_package.py --privacy
```

`--privacy` requires a nonempty token list, `pypdf`, Pillow, and Tesseract. It adds
encoded-token checks, PDF object/text/metadata inspection, and image metadata/OCR
checks. It does not OCR pixels embedded in PDFs; image-only pages without
extracted text are rejected for separate review. An ordinary manifest pass is
not a substitute for this privacy scan, and no scan guarantees unlinkability.

To create the local repository from an unchanged extracted package:

```sh
python3 tools/initialize_anonymous_git.py
```

The initializer disables signing and hooks for its initial commit and verifies
the committed bytes. It does not configure a remote, contact GitHub, or publish
anything. GitHub account ownership remains separate from commit identity; follow
[PUBLICATION.md](../PUBLICATION.md) before publishing.

The saved release reports are [PACKAGE_CHECKS.json](../PACKAGE_CHECKS.json) for
bounded scientific checks and [PRIVACY_REPORT.json](../PRIVACY_REPORT.json) for
the privacy audit. Adding or changing files makes the original manifest check
fail; preserve the delivered package as your baseline.

After editing or adding research files, maintain the navigation with:

```sh
python3 tools/build_navigation.py
```

This updates `docs/catalog/` and `docs/REFERENCES.md`. It does not update the
integrity manifest or certify a changed release. Review the new navigation and
repeat the release checks before distributing a modified package.
