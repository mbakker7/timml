import os
import tempfile
from pathlib import Path

import nbformat
import pytest
from nbconvert.preprocessors import ExecutePreprocessor

nbdirs = [
    os.path.join("docs/00userguide"),
    os.path.join("docs/02examples"),
    os.path.join("docs/03xsections"),
    os.path.join("docs/04tests"),
]

testdir = tempfile.mkdtemp()


def get_notebooks():
    skip = ["benchmarking_besselaes.ipynb", "anisotropy_pap_posted.ipynb"]
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


# def get_jupyter_kernel():
#     try:
#         jklcmd = ("jupyter", "kernelspec", "list")
#         b = subprocess.Popen(
#             jklcmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
#         ).communicate()[0]
#         if isinstance(b, bytes):
#             b = b.decode("utf-8")
#         print(b)
#         for line in b.splitlines():
#             if "python" in line:
#                 kernel = line.split()[0]
#     except:  # noqa
#         kernel = None

#     return kernel


# @pytest.mark.notebooks
# @pytest.mark.parametrize("pth", get_notebooks())
# def test_notebook(pth):
#     kernel = get_jupyter_kernel()
#     print("available jupyter kernel {}".format(kernel))

#     fn = os.path.basename(pth)

#     cmd = (
#         "jupyter "
#         + "nbconvert "
#         + "--ExecutePreprocessor.timeout=600 "
#         + "--to "
#         + "notebook "
#         + "--execute "
#         + "{} ".format(pth)
#         + "--output-dir "
#         + "{} ".format(testdir)
#         + "--output "
#         + "{}".format(fn)
#     )
#     ival = os.system(cmd)
#     assert ival == 0, "could not run {}".format(fn)


# %%
if __name__ == "__main__":
    for notebook in get_notebooks():
        os.system(
            "jupyter nbconvert --clear-output --inplace "
            f"--ClearMetadataPreprocessor.enabled=True {notebook}"
        )
