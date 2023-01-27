import pytest

def pytest_addoption(parser):
    parser.addoption("--my_test_dump", action="store", default=None,
        help="Print test items in my custom format")

def pytest_collection_finish(session):
    if session.config.option.my_test_dump is not None:
        for item in session.items:
            print('{}::{}'.format(item.fspath, item.name))
        pytest.exit('Done!')
    # make this do the cargo tests and comparison too, and not keep files around
