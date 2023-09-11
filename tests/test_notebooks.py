import os
import shutil
import subprocess
import tempfile

import pytest

nbdir = os.path.join("notebooks")

testdir = tempfile.mkdtemp()


def get_notebooks():
    skip = ["benchmarking_besselaes.ipynb"]
    return [f for f in os.listdir(nbdir) if f.endswith(".ipynb") and f not in skip]


def get_jupyter_kernel():
    try:
        jklcmd = ("jupyter", "kernelspec", "list")
        b = subprocess.Popen(
            jklcmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        ).communicate()[0]
        if isinstance(b, bytes):
            b = b.decode("utf-8")
        print(b)
        for line in b.splitlines():
            if "python" in line:
                kernel = line.split()[0]
    except:
        kernel = None

    return kernel


@pytest.mark.notebooks
@pytest.mark.parametrize("fn", get_notebooks())
def test_notebook(fn):
    kernel = get_jupyter_kernel()
    print("available jupyter kernel {}".format(kernel))

    pth = os.path.join(nbdir, fn)

    cmd = (
        "jupyter "
        + "nbconvert "
        + "--ExecutePreprocessor.timeout=600 "
        + "--to "
        + "notebook "
        + "--execute "
        + "{} ".format(pth)
        + "--output-dir "
        + "{} ".format(testdir)
        + "--output "
        + "{}".format(fn)
    )
    ival = os.system(cmd)
    assert ival == 0, "could not run {}".format(fn)


if __name__ == "__main__":
    test_notebook()
    shutil.rmtree(testdir)
