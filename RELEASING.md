# Releasing TCRsift

This document explains how to publish a new TCRsift release once your changes are merged and the working tree is clean.

## Fast path (recommended)

Use the release script to run lint/tests, optionally bump the version, build, upload, and tag:

```bash
# Verify readiness without uploading or tagging
./deploy.sh --dry-run

# Bump to a new version and release
./deploy.sh 0.2.5

# Release the current version without bumping
./deploy.sh
```

## Notes

- Version is stored in `tcrsift/version.py` as `__version__`.
- The release script will refuse to run if the working tree is dirty.
