# Releasing

## Build Locally

```bash
uv sync --extra dev
uv build
```

Artifacts are written to `dist/`.

## CI

The repository includes:

- a CI workflow that installs dependencies, runs tests, and builds the package
- a release workflow that builds distributions on tags and can publish to PyPI if credentials are configured

## Suggested Release Process

1. Update documentation and bundled data notes if needed.
2. Run `uv run pytest`.
3. Run `uv build`.
4. Create and push a version tag such as `v0.9.1`.
5. Let GitHub Actions build artifacts and publish if PyPI publishing is configured.
