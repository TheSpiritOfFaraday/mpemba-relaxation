#!/usr/bin/env python3
"""Verify release bytes and, optionally, private identity and media hygiene.

This verifies the package as delivered, before running programs that regenerate
its receipts or plots. It never writes files. A successful scan cannot establish
that distinctive research or a publishing account is untraceable.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
from urllib.parse import quote

sys.dont_write_bytecode = True
MANIFEST = "PACKAGE_MANIFEST.json"
MAX_FILE_BYTES = 100 * 1024 * 1024
IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".webp", ".tif", ".tiff", ".bmp"}
FORBIDDEN_PARTS = {".DS_Store", "__MACOSX", "__pycache__", ".venv", "node_modules", ".lake", ".tools"}
FORBIDDEN_NAMES = {".env", ".netrc", ".npmrc", ".pypirc", "credentials", "credentials.json", "id_rsa", "id_ed25519"}
FORBIDDEN_SUFFIXES = {".pyc", ".pyo", ".p12", ".pfx", ".pem", ".key"}
# Construct path expressions in pieces so the scanner does not flag its source.
HOME_PATTERN = re.compile(
    r"/(?:" + "Users" + "|" + "home" + r")/(?!anonymous(?:[-/]|\b))[^\s/<>\"'\\]+",
    re.I,
)
SECRET_PATTERNS = [
    re.compile(r"gh[pousr]_[A-Za-z0-9]{20,}"),
    re.compile(r"github_pat_[A-Za-z0-9_]{30,}"),
    re.compile(r"AKIA[0-9A-Z]{16}"),
    re.compile(r"xox[baprs]-[A-Za-z0-9-]{20,}"),
    re.compile(r"sk-(?:proj-)?[A-Za-z0-9_-]{32,}"),
    re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----"),
]


class VerificationError(Exception):
    """A controlled error whose message contains no scanned secret values."""


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise VerificationError("JSON contains duplicate object keys")
        result[key] = value
    return result


def load_tokens():
    try:
        tokens = json.loads(os.environ.get("ANON_PRIVATE_TOKENS", "[]"))
    except (ValueError, TypeError):
        raise VerificationError("ANON_PRIVATE_TOKENS must be a JSON list of private strings") from None
    if not isinstance(tokens, list) or any(not isinstance(t, str) or len(t) < 3 for t in tokens):
        raise VerificationError("Private token entries must be strings of at least three characters")
    return list(dict.fromkeys(tokens))


def safe_path(path, tokens):
    return "[redacted path]" if any(t.casefold() in path.casefold() for t in tokens) else path


def variants(tokens):
    result = set()
    for token in tokens:
        for encoding in ("utf-8", "utf-16le", "utf-16be"):
            result.add(token.encode(encoding).lower())
        result.add(base64.b64encode(token.encode()).lower())
        result.add(quote(token, safe="").encode().lower())
    return result


def inspect_text(text, label, tokens):
    if any(t.casefold() in text.casefold() for t in tokens):
        raise VerificationError(f"Private identifier detected: {label}")
    if HOME_PATTERN.search(text):
        raise VerificationError(f"Personal home-directory path detected: {label}")
    if any(pattern.search(text) for pattern in SECRET_PATTERNS):
        raise VerificationError(f"Credential-like content detected: {label}")


def inspect_bytes(data, label, tokens, encoded_tokens):
    if any(token in data.lower() for token in encoded_tokens):
        raise VerificationError(f"Encoded private identifier detected: {label}")
    # Binary objects are also examined without assuming valid UTF-8.
    inspect_text(data.decode("utf-8", "ignore"), label, tokens)


def pdf_object_strings(reader):
    seen = set()
    stack = [reader.trailer]
    while stack:
        obj = stack.pop()
        if hasattr(obj, "get_object"):
            obj = obj.get_object()
        marker = id(obj)
        if marker in seen:
            continue
        seen.add(marker)
        if isinstance(obj, str):
            yield obj.encode("utf-8")
        elif isinstance(obj, bytes):
            yield obj
        elif isinstance(obj, dict):
            stack.extend(obj.keys())
            stack.extend(obj.values())
            if hasattr(obj, "get_data"):
                yield obj.get_data()
        elif isinstance(obj, (list, tuple)):
            stack.extend(obj)


def inspect_pdf(data, label, tokens, encoded_tokens):
    from pypdf import PdfReader
    reader = PdfReader(io.BytesIO(data))
    if reader.is_encrypted:
        raise VerificationError(f"Encrypted PDF cannot be audited: {label}")
    if reader.metadata or "/Metadata" in reader.trailer["/Root"] or "/ID" in reader.trailer:
        raise VerificationError(f"PDF document metadata remains: {label}")
    if reader.attachments:
        raise VerificationError(f"PDF contains embedded attachments: {label}")
    for raw in pdf_object_strings(reader):
        inspect_bytes(raw, label + " (PDF object)", tokens, encoded_tokens)
    for page in reader.pages:
        extracted = page.extract_text() or ""
        if not extracted.strip() and len(page.images):
            raise VerificationError(f"PDF image-only page requires separate OCR review: {label}")
        inspect_text(extracted, label + " (PDF text)", tokens)
    return len(reader.pages)


def inspect_image(path, label, tokens, encoded_tokens, tesseract):
    from PIL import Image
    with Image.open(path) as image:
        inspect_text(repr(image.info) + repr(dict(image.getexif())), label + " (image metadata)", tokens)
        # Scan every frame when animation or multiple image pages are present.
        frame_count = getattr(image, "n_frames", 1)
        for index in range(frame_count):
            image.seek(index)
            buffer = io.BytesIO()
            image.convert("RGB").save(buffer, format="PNG")
            run = subprocess.run(
                [tesseract, "stdin", "stdout"], input=buffer.getvalue(),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120,
                check=False,
            )
            if run.returncode:
                raise VerificationError(f"Image OCR did not complete: {label}")
            inspect_bytes(run.stdout, label + " (image OCR)", tokens, encoded_tokens)
    return frame_count


def validate_relative_path(rel):
    if not isinstance(rel, str) or not rel or "\\" in rel or any(ord(c) < 32 for c in rel):
        raise VerificationError("Manifest contains an unsafe relative path")
    parts = PurePosixPath(rel).parts
    if rel.startswith("/") or any(p in {"", ".", "..", ".git"} for p in rel.split("/")):
        raise VerificationError("Manifest contains an unsafe relative path")
    if PurePosixPath(rel).as_posix() != rel:
        raise VerificationError("Manifest path is not normalized")
    if any(p in FORBIDDEN_PARTS or p in FORBIDDEN_NAMES or p.startswith("._") for p in parts):
        raise VerificationError("Manifest contains a private, cache, or filesystem-metadata file")
    if PurePosixPath(rel).suffix.lower() in FORBIDDEN_SUFFIXES or PurePosixPath(rel).name.startswith(".env."):
        raise VerificationError("Manifest contains a private or compiled file")


def verify(root, privacy=False):
    root = Path(root).resolve()
    tokens = load_tokens()
    encoded_tokens = variants(tokens)
    if privacy and not tokens:
        raise VerificationError("--privacy requires ANON_PRIVATE_TOKENS with the identifiers to check")
    tesseract = None
    if privacy:
        try:
            import pypdf  # noqa: F401
            import PIL  # noqa: F401
        except ImportError:
            raise VerificationError("--privacy requires pypdf and Pillow") from None
        tesseract = shutil.which("tesseract")
        if not tesseract:
            raise VerificationError("--privacy requires tesseract for image OCR")
    manifest_path = root / MANIFEST
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise VerificationError("A regular PACKAGE_MANIFEST.json is required")
    manifest_bytes = manifest_path.read_bytes()
    inspect_bytes(manifest_bytes, MANIFEST, tokens, encoded_tokens)
    try:
        manifest = json.loads(manifest_bytes, object_pairs_hook=unique_object)
    except (ValueError, UnicodeDecodeError):
        raise VerificationError("The package manifest is not valid JSON") from None
    inventory = manifest.get("files") if isinstance(manifest, dict) else None
    if not isinstance(inventory, dict) or not inventory:
        raise VerificationError("Manifest files must be a nonempty path-to-SHA256 object")
    if MANIFEST in inventory:
        raise VerificationError("The manifest must not include its own hash")
    for rel, digest in inventory.items():
        validate_relative_path(rel)
        if not isinstance(digest, str) or not re.fullmatch(r"[a-f0-9]{64}", digest):
            raise VerificationError("Manifest contains an invalid SHA256 digest")
    actual = set()
    for parent, dirs, names in os.walk(root, followlinks=False):
        for directory in dirs[:]:
            path = Path(parent) / directory
            if path.is_symlink():
                raise VerificationError("Package contains a directory symlink")
            if path == root / ".git":
                dirs.remove(directory)
        for name in names:
            path = Path(parent) / name
            if path.is_symlink() or not path.is_file():
                raise VerificationError("Package contains a symlink or nonregular file")
            rel = path.relative_to(root).as_posix()
            if rel != MANIFEST:
                validate_relative_path(rel)
                actual.add(rel)
    if actual != set(inventory):
        raise VerificationError(
            "Package file inventory differs from manifest "
            f"({len(set(inventory) - actual)} missing; {len(actual - set(inventory))} extra)"
        )
    counts = {"files": 0, "bytes": 0, "pdfs": 0, "pdf_pages": 0, "images": 0, "image_frames": 0}
    for rel, expected in sorted(inventory.items()):
        path = root / rel
        label = safe_path(rel, tokens)
        inspect_text(rel, label + " (filename)", tokens)
        if path.stat().st_size > MAX_FILE_BYTES:
            raise VerificationError(f"File exceeds the 100 MiB GitHub limit: {label}")
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            raise VerificationError(f"SHA256 mismatch: {label}")
        inspect_bytes(data, label, tokens, encoded_tokens)
        counts["files"] += 1
        counts["bytes"] += len(data)
        if privacy and path.suffix.lower() == ".pdf":
            counts["pdf_pages"] += inspect_pdf(data, label, tokens, encoded_tokens)
            counts["pdfs"] += 1
        elif privacy and path.suffix.lower() in IMAGE_SUFFIXES:
            counts["image_frames"] += inspect_image(path, label, tokens, encoded_tokens, tesseract)
            counts["images"] += 1
        elif privacy:
            try:
                data.decode("utf-8")
            except UnicodeDecodeError:
                raise VerificationError(f"Unsupported binary file requires a separate audit: {label}") from None
    return {
        "passed": True,
        "scope": "manifest, file inventory, SHA256, paths, known credential patterns, size limits"
        + (", supplied private identifiers, PDF objects/text/metadata, image metadata and OCR" if privacy else ""),
        "counts": counts,
        "private_identifier_count": len(tokens),
        "privacy_media_checks_run": privacy,
        "limits": "Does not guarantee anonymity of distinctive content, prior disclosures, account activity, or later edits. PDF image pixels are not OCR-scanned; image-only pages without extracted text are rejected.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--privacy", action="store_true")
    args = parser.parse_args()
    try:
        report = verify(args.root, args.privacy)
    except VerificationError as error:
        print(json.dumps({"passed": False, "error": str(error)}, indent=2))
        return 1
    except Exception as error:
        # Uncontrolled library exceptions can contain private paths or strings.
        print(json.dumps({"passed": False, "error": "Audit failed with " + type(error).__name__}, indent=2))
        return 1
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
