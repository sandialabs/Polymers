import pytest
from re import sub
from subprocess import run


def pytest_addoption(parser):
    parser.addoption(
        "--compare", action="store", default=None,
        help="Print test items in my custom format"
    )


def pytest_collection_finish(session):
    f = open("__pycache__/cargo.tests", "w")
    run(
        ['cargo', 'test', '--tests', '--', '--list', '--format=terse'],
        stdout=f, stderr=0
    )
    f.close()
    run(
        ['sed', '-i', 's@: test@@', '__pycache__/cargo.tests']
    )
    f = open("__pycache__/julia.tests", "w")
    run(
        ['grep', '-r', '@testset', 'src/'],
        stdout=f, stderr=0
    )
    f.close()
    run(
        ['sed', '-i', 's@^.*testset "@@', '__pycache__/julia.tests']
    )
    run(
        ['sed', '-i', 's@".*$@@', '__pycache__/julia.tests']
    )
    stdout = open('__pycache__/pytest.tests', 'w')
    if session.config.option.compare is not None:
        for item in session.items:
            class_name = sub(
                r'(?<!^)(?=[A-Z])', '_', item.parent.name
            ).lower()
            stdout.write(
                '{}::{}::{}\n'.format(
                    item.fspath, class_name, item.name
                )
            )
        stdout.close()
        run(
            ['sed', '-i', 's@.py::@::@', '__pycache__/pytest.tests']
        )
        run(
            ['sed', '-i', 's@test_@@', '__pycache__/pytest.tests']
        )
        run(
            ['sed', '-i', 's@^.*src/@@', '__pycache__/pytest.tests']
        )
        run(
            ['sed', '-i', 's@/@::@g', '__pycache__/pytest.tests']
        )
        run(
            ['sort', '__pycache__/cargo.tests',
             '-o', '__pycache__/cargo.tests']
        )
        run(
            ['sort', '__pycache__/pytest.tests',
             '-o', '__pycache__/pytest.tests']
        )
        run(
            ['sort', '__pycache__/julia.tests',
             '-o', '__pycache__/julia.tests']
        )
        code = run(
            ['cmp', '-s', '__pycache__/cargo.tests',
             '__pycache__/pytest.tests']
        ).returncode
        code += run(
            ['cmp', '-s', '__pycache__/cargo.tests',
             '__pycache__/julia.tests']
        ).returncode
        if code == 0:
            pytest.exit('tests match across languages', code)
        else:
            pytest.exit('tests do not match across languages', code)
