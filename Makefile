.PHONY: install test check

install:
	uv sync --extra notebooks --extra test

test:
	uv run --extra test python -m pytest

check: test
	uv run python -m compileall -q scripts tests
