from __future__ import division, absolute_import, print_function

import sys

def test_import():
    try:
        import timml
    except:
        fail = True
        assert fail is False, 'could not import timml'
    return

if __name__ == '__main__':
    test_import()
