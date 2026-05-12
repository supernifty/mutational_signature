# Development

## Recommended Workflow

```bash
uv sync --extra dev
uv run pytest
uv build
```

`uv` is the preferred contributor workflow because it gives fast environment creation and repeatable dependency resolution. The package still supports standard `pip` installation for users who do not want to adopt `uv`.

## Compatibility Workflow

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .[dev]
pytest
python -m build
```

## Repository Standards

- Keep changes targeted. Many modules are exposed directly as CLI tools.
- Match the existing Python style in touched files. Several modules use 2-space indentation.
- Prefer adding small synthetic fixtures over relying on large external datasets.
- Document any changes to bundled reference signature files, including source and version.

## Validation

At minimum, verify:

- `count` on the small example fixture
- `decompose` on the resulting count file
- `plot-counts` on the resulting count file

The automated smoke tests cover that exact path.
