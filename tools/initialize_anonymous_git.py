#!/usr/bin/env python3
"""Create one locally verified anonymous commit; never configure or contact a remote.

Run only in an unchanged extracted release. Existing Git repositories are refused.
Publishing still requires a separately chosen account with no identifying links.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
from verify_package import VerificationError, verify  # noqa: E402

AUTHOR = "Anonymous Researcher"
EMAIL = "anonymous@example.invalid"
DATE = "2000-01-01T00:00:00+0000"
MESSAGE = "Initial anonymous research release"


def main():
    root = Path(__file__).resolve().parents[1]
    if os.path.lexists(root / ".git"):
        raise VerificationError("Refusing to alter an existing Git repository")
    verify(root)
    git = shutil.which("git")
    if not git:
        raise VerificationError("Git is required")
    # Remove inherited Git controls, author values, signing settings and
    # configuration injection. Explicit config paths disable personal config.
    env = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    with tempfile.TemporaryDirectory(prefix="anonymous-git-") as temporary:
        neutral_directory = Path(temporary)
        template = neutral_directory / "empty-template"
        template.mkdir()
        env.update({
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_CONFIG_GLOBAL": os.devnull,
            "GIT_CONFIG_SYSTEM": os.devnull,
            "GIT_ATTR_NOSYSTEM": "1",
            "GIT_TERMINAL_PROMPT": "0",
            "GIT_AUTHOR_NAME": AUTHOR,
            "GIT_AUTHOR_EMAIL": EMAIL,
            "GIT_COMMITTER_NAME": AUTHOR,
            "GIT_COMMITTER_EMAIL": EMAIL,
            "GIT_AUTHOR_DATE": DATE,
            "GIT_COMMITTER_DATE": DATE,
            "TZ": "UTC",
            "LC_ALL": "C",
        })

        def run(*arguments):
            result = subprocess.run([git, *arguments], cwd=root, env=env, capture_output=True, text=True)
            if result.returncode:
                raise VerificationError("Git operation failed; inspect the new local repository before retrying")
            return result.stdout.strip()

        run("init", "--initial-branch=main", "--template=" + str(template), ".")
        settings = {
            "user.name": AUTHOR,
            "user.email": EMAIL,
            "user.useConfigOnly": "true",
            "commit.gpgsign": "false",
            "tag.gpgsign": "false",
            "core.logAllRefUpdates": "false",
            "core.hooksPath": ".git/disabled-hooks",
            "core.autocrlf": "false",
            "core.attributesFile": os.devnull,
            "core.excludesFile": os.devnull,
        }
        for key, value in settings.items():
            run("config", "--local", key, value)
        (root / ".git" / "disabled-hooks").mkdir()
        # Force-add precisely the manifest, unaffected by .gitignore rules.
        manifest = json.loads((root / "PACKAGE_MANIFEST.json").read_text())
        names = sorted(manifest["files"]) + ["PACKAGE_MANIFEST.json"]
        for start in range(0, len(names), 200):
            run("add", "--force", "--", *names[start:start + 200])
        run("-c", "core.hooksPath=.git/disabled-hooks", "commit", "--no-gpg-sign", "-m", MESSAGE)
        author = run("show", "-s", "--format=%an%n%ae%n%cn%n%ce%n%ai%n%ci", "HEAD").splitlines()
        expected = [AUTHOR, EMAIL, AUTHOR, EMAIL, "2000-01-01 00:00:00 +0000", "2000-01-01 00:00:00 +0000"]
        if author != expected:
            raise VerificationError("Commit identity/date verification failed")
        if run("rev-list", "--count", "HEAD") != "1" or run("remote"):
            raise VerificationError("Repository history or remote verification failed")
        tree = run("ls-tree", "-r", "--name-only", "HEAD").splitlines()
        if set(tree) != set(names):
            raise VerificationError("Committed file inventory differs from the verified package")
        committed_objects = {}
        for entry in run("ls-tree", "-r", "-z", "HEAD").split("\0"):
            if entry:
                descriptor, name = entry.split("\t", 1)
                mode, kind, object_id = descriptor.split()
                if kind != "blob" or mode not in {"100644", "100755"}:
                    raise VerificationError("Unexpected object type in committed package")
                committed_objects[name] = object_id
        for start in range(0, len(names), 200):
            batch = names[start:start + 200]
            hashes = run("hash-object", "--no-filters", "--", *batch).splitlines()
            if len(hashes) != len(batch) or any(committed_objects[name] != digest for name, digest in zip(batch, hashes)):
                raise VerificationError("Git transformed package bytes during staging")
        if "gpgsig" in run("cat-file", "-p", "HEAD"):
            raise VerificationError("Unexpected commit signature")
        if run("status", "--porcelain", "--untracked-files=all"):
            raise VerificationError("Repository is not clean after initialization")
        verify(root)
        print(json.dumps({
            "passed": True,
            "commit": run("rev-parse", "HEAD"),
            "author": AUTHOR,
            "email": EMAIL,
            "commit_date": "2000-01-01T00:00:00Z",
            "branch": "main",
            "committed_files": len(names),
            "remotes": [],
            "published": False,
        }, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except VerificationError as error:
        print(json.dumps({"passed": False, "error": str(error)}, indent=2))
        raise SystemExit(1)
    except Exception as error:
        print(json.dumps({"passed": False, "error": "Initialization failed with " + type(error).__name__}, indent=2))
        raise SystemExit(1)
