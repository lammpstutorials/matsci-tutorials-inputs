#!/usr/bin/env python3
"""Validate the tutorial file trees against their .manifest files.

Works for both tutorial repository layouts:
- the article repository:  files/tutorial<N>/ with .manifest and solution/
- the download repository: tutorial<N>/ with .manifest and solution/

The script enforces what the LAMMPS-GUI tutorial wizard actually does when it
downloads a tutorial (see setupTutorial() / downloadTutorialFiles() in
lammps-gui/src/lammpsgui.cpp):

- every .manifest entry is fetched from the repository, so a listed file that
  does not exist produces an HTTP 404 error and aborts the download (ERROR);
- files that exist but are not listed are never downloaded (WARNING);
- only entries in the solution/ subfolder are downloaded from subfolders, and
  only when the user requests the solutions; entries in any other subfolder
  are silently ignored by the GUI (ERROR);
- the first top-level entry is opened in the editor after the download, so a
  manifest without any top-level entry downloads but opens nothing (ERROR);
- a symbolic link is downloaded as a small placeholder file and the GUI then
  replaces it with a copy of the file it points to.  This only works when the
  link points to "../<same file name>" and that file is listed in the
  manifest *before* the link, so it has already been downloaded (ERROR
  otherwise).  A top-level symlink would point outside the download
  directory and can never be resolved (ERROR).

Exit status: 0 when no errors were found (warnings do not fail the check),
1 when at least one error was found, 2 on usage problems.
"""

import argparse
import re
import sys
from pathlib import Path

TUTORIAL_DIR = re.compile(r"^tutorial(\d+)$")


def parse_manifest(path):
    """Return the manifest entries with their line numbers, skipping comments."""
    entries = []
    for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        entries.append((lineno, line))
    return entries


def check_tutorial(tdir):
    """Validate one tutorial directory; return (errors, warnings) message lists."""
    errors = []
    warnings = []

    manifest = tdir / ".manifest"
    if not manifest.is_file():
        return (["no .manifest file"], [])

    entries = parse_manifest(manifest)
    seen = {}
    toplevel_before = set()  # top-level entries listed so far, for symlink order
    have_toplevel = False

    for lineno, entry in entries:
        where = f".manifest:{lineno}: {entry}"

        if entry.startswith("/") or "\\" in entry or ".." in entry.split("/"):
            errors.append(f"{where}: invalid path (absolute, backslash, or '..')")
            continue
        if entry in seen:
            errors.append(f"{where}: duplicate entry (first listed on line {seen[entry]})")
            continue
        seen[entry] = lineno

        parts = entry.split("/")
        if len(parts) > 1 and parts[0] != "solution":
            errors.append(f"{where}: subfolder entries other than solution/ are "
                          "ignored by the GUI and never downloaded")
            continue
        if len(parts) > 2:
            warnings.append(f"{where}: nested more than one level below solution/")

        fpath = tdir / entry
        if not fpath.exists() and not fpath.is_symlink():
            errors.append(f"{where}: listed but missing on disk "
                          "(download would fail with HTTP 404)")
            continue

        if fpath.is_symlink():
            target = fpath.readlink().as_posix()
            if len(parts) == 1:
                errors.append(f"{where}: top-level symlink; its placeholder points "
                              "outside the download directory and is never resolved")
            elif target != f"../{fpath.name}":
                errors.append(f"{where}: symlink target '{target}' is not "
                              f"'../{fpath.name}'; the GUI only resolves links to the "
                              "same file name in the tutorial folder")
            elif not fpath.resolve().exists():
                errors.append(f"{where}: dangling symlink (target does not exist)")
            elif fpath.name not in toplevel_before:
                errors.append(f"{where}: link target '{fpath.name}' must be listed "
                              "before this entry so it is downloaded first")
        if len(parts) == 1:
            have_toplevel = True
            toplevel_before.add(entry)

    if entries and not have_toplevel:
        errors.append(".manifest: no top-level entry; the GUI would open nothing "
                      "after the download")
    if not entries:
        errors.append(".manifest: empty")

    # files on disk that the manifest does not cover
    for fpath in sorted(tdir.rglob("*")):
        if not fpath.is_file() and not fpath.is_symlink():
            continue
        rel = fpath.relative_to(tdir).as_posix()
        if rel == ".manifest" or fpath.name.startswith("."):
            continue
        if rel not in seen:
            warnings.append(f"{rel}: on disk but not in .manifest "
                            "(will not be downloaded)")

    return (errors, warnings)


def find_tutorials(root):
    """Return the tutorial directories under root, handling both layouts."""
    base = root / "files" if (root / "files").is_dir() else root
    found = []
    for path in base.iterdir():
        m = TUTORIAL_DIR.match(path.name)
        if path.is_dir() and m:
            found.append((int(m.group(1)), path))
    return base, sorted(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("root", nargs="?", default=".", type=Path,
                        help="repository root to check (default: current directory)")
    args = parser.parse_args()

    if not args.root.is_dir():
        parser.error(f"not a directory: {args.root}")
    base, tutorials = find_tutorials(args.root)
    if not tutorials:
        parser.error(f"no tutorial<N> directories found under {base}")

    numbers = [num for num, _ in tutorials]
    total_errors = total_warnings = 0

    for num, tdir in tutorials:
        errors, warnings = check_tutorial(tdir)
        total_errors += len(errors)
        total_warnings += len(warnings)
        if errors or warnings:
            print(f"{tdir.relative_to(args.root)}:")
            for msg in errors:
                print(f"  ERROR:   {msg}")
            for msg in warnings:
                print(f"  WARNING: {msg}")

    missing = sorted(set(range(1, max(numbers) + 1)) - set(numbers))
    if missing:
        print(f"note: missing tutorial numbers: {', '.join(map(str, missing))}")

    print(f"checked {len(tutorials)} tutorials: "
          f"{total_errors} errors, {total_warnings} warnings")
    return 1 if total_errors else 0


if __name__ == "__main__":
    sys.exit(main())
