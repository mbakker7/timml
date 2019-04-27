def test_import():
    try:
        import timml
    except:
        fail = True
        assert fail is False, "could not import timml"
    return


if __name__ == "__main__":
    test_import()
