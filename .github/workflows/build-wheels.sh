#!/bin/bash
set -e -x

for PYBIN in /opt/python/cp3[7891]*/bin; do
    "${PYBIN}/pip" install cffi maturin
    "${PYBIN}/maturin" build --features python -i "${PYBIN}/python" --release
done

for wheel in target/wheels/*.whl; do
    auditwheel repair "${wheel}"
done
