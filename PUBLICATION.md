# Publishing anonymously

This package removes identifying content and metadata within the checked scope.
It cannot guarantee anonymity. A GitHub account, its public activity, reused
identifiers, previously disclosed copies, and distinctive research content can
link a publication to a person. GitHub also holds account and access information.

The ZIP contains no Git history. Do not publish the original working directory,
an existing repository, a local backup, or an account-linked conversation.

## Prepare the publishing identity

Use a separate GitHub account with a nonidentifying username, profile, avatar,
and contact information. Review that account's visible activity and account
settings before publication. GitHub email-privacy settings alone do not make a
known account anonymous. A GitHub `noreply` address may still contain an account
identifier; the package's Git helper uses a neutral author identity instead.

Use a clean extraction of this package. Inspect the privacy report and files,
verify the published inventory, then initialize a fresh repository with:

```sh
python3 tools/verify_package.py
python3 tools/initialize_anonymous_git.py
git status
```

The helper creates one initial commit on `main`, sets a fixed anonymous author,
email, and neutral commit date, and disables Git signing and hooks for this
repository. It refuses an existing Git repository and configures no remote.
Follow its output to inspect the initialized state.
Do not replace the neutral identity with a personal email, import old commits,
enable signing with a personal key, or reuse an existing remote.

Before publishing any commit, inspect its author and committer:

```sh
git log --format=fuller
git remote -v
git config --local --get user.name
git config --local --get user.email
```

Commit metadata and repository ownership are separate: an anonymous author does
not anonymize the GitHub account that owns or pushes the repository.

## Deliberate publication

The package preparation itself does not create or publish a GitHub repository.
Run these steps only after checking the exact files and selecting the intended
separate account. With GitHub CLI installed, inspect its active account:

```sh
gh auth status --hostname github.com
```

If it is the wrong account, sign in to the intended one before proceeding. Replace
`ANONYMOUS-ACCOUNT` below with that account's nonidentifying username. Create the
repository without copying files from any other workspace:

```sh
gh repo create ANONYMOUS-ACCOUNT/mpemba-relaxation --public
git remote add origin https://github.com/ANONYMOUS-ACCOUNT/mpemba-relaxation.git
git push --set-upstream origin main
```

Inspect the public repository, file contents, commit authors/committers, profile,
and rendered PDFs afterward. Avoid public links to personal profiles, accounts,
cloud storage, old releases, or conversations.

## Later updates

Re-run the package privacy verification before every public update. Reproduction
can write local paths into receipts and author or tool metadata into PDFs. New
archives, images, datasets, or logs need their own checks. Do not assume that an
earlier successful scan covers later generated files.

The optional `python3 tools/verify_package.py --privacy` mode also checks PDF
objects/text/metadata and image metadata/OCR. It requires the document dependencies
and Tesseract. Private identifiers can be supplied through the
`ANON_PRIVATE_TOKENS` JSON-list environment variable; never save the private list
in the repository or a public command transcript. The ordinary check validates
the release manifest and known credential patterns; it is not a substitute for
media checks or a guarantee about personal details unknown to the scanner.

No license has been selected automatically. Review [RIGHTS.md](RIGHTS.md) before
adding a license or claiming open-source status.

GitHub's documentation explains [commit email attribution](https://docs.github.com/en/account-and-profile/how-tos/email-preferences/setting-your-commit-email-address)
and the [limits of removing previously published sensitive data](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/removing-sensitive-data-from-a-repository).
