import os
import tempfile
from pathlib import Path

import nbformat
import pytest
from nbconvert.preprocessors import (
    ClearMetadataPreprocessor,
    ClearOutputPreprocessor,
    ExecutePreprocessor,
)

nbdirs = [
    os.path.join("docs/00userguide/tutorials"),
    os.path.join("docs/00userguide/howtos"),
    os.path.join("docs/02examples"),
    os.path.join("docs/03xsections"),
    os.path.join("docs/04tests"),
]

testdir = tempfile.mkdtemp()


def get_notebooks():
    skip = ["benchmarking_besselaes.ipynb", "vertical_anisotropy.ipynb"]
    nblist = []
    for nbdir in nbdirs:
        nblist += [
            os.path.join(nbdir, f)
            for f in os.listdir(nbdir)
            if f.endswith(".ipynb") and f not in skip
        ]
    # nblist = ["notebooks/test_plot.ipynb"]
    return nblist


# @pytest.mark.notebooks
@pytest.mark.skip(reason="Use pytest --nbval on notebooks directly for coverage.")
@pytest.mark.parametrize("pth", get_notebooks())
def test_notebook_py(pth):
    pth = Path(pth)
    with open(pth, "r", encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)
        ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
        try:
            assert ep.preprocess(nb, {"metadata": {"path": pth.parent}}) is not None, (
                f"Got empty notebook for {os.path.basename(pth)}"
            )
        except Exception as e:
            pytest.fail(reason=f"Failed executing {os.path.basename(pth)}: {e}")


# %% clear output and metadata of all notebooks
if __name__ == "__main__":
    clear_output = ClearOutputPreprocessor()
    clear_metadata = ClearMetadataPreprocessor()

    for notebook in get_notebooks():
        print("Clearing notebook:", notebook)
        with open(notebook, "r", encoding="utf-8") as f:
            nb = nbformat.read(f, as_version=4)

        # run nbconvert preprocessors to clear outputs and metadata
        clear_output.preprocess(nb, {})
        clear_metadata.preprocess(nb, {})

        with open(notebook, "w", encoding="utf-8") as f:
            nbformat.write(nb, f)
