import pytest
from sys import stdout
from subprocess import run

def pytest_addoption(parser):
    parser.addoption("--compare", action="store", default=None,
        help="Print test items in my custom format")

def pytest_collection_finish(session):
    f = open("__pycache__/cargo.tests", "w")
    run(['cargo', 'test', '--tests', '--', '--list', '--format=terse'], stdout=f)
    f.close()
    stdout = open('__pycache__/pytest.tests', 'w')
    if session.config.option.compare is not None:
        for item in session.items:
            stdout.write('{}::{}\n'.format(item.fspath, item.name))
        stdout.close()
        code = run(['cmp', '__pycache__/cargo.tests', '__pycache__/pytest.tests']).returncode
        if code == 0:
            pytest.exit('tests match across languages', code)
        else:
            pytest.exit('tests do not match across languages', code)
