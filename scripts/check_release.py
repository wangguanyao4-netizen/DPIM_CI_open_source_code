"""Static release checks only: no imports of experiment or plotting modules."""
from pathlib import Path
import ast
import hashlib
import json
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'provenance/release_manifest.json'


def main():
    records = json.loads(MANIFEST.read_text(encoding='utf-8'))['files']
    checked_python = 0
    excluded = {'.mat', '.npz', '.npy', '.h5', '.hdf5', '.pdf', '.png', '.jpg', '.zip'}
    data_manifest = json.loads((ROOT / 'provenance/data_sources.json').read_text(encoding='utf-8'))
    approved_data = {item['path'] for item in data_manifest['files']}
    actual_data = {p.relative_to(ROOT).as_posix() for p in (ROOT / 'data').rglob('*') if p.is_file() and p.suffix != '.md'}
    assert actual_data == approved_data
    for item in records:
        path = (ROOT / item['path']).resolve()
        assert path.is_relative_to(ROOT.resolve())
        assert path.is_file(), item['path']
        assert path.suffix.lower() not in excluded, item['path']
        if path.suffix.lower() == '.csv':
            assert item['path'] in approved_data, item['path']
        raw = path.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == item['sha256'], item['path']
        if path.suffix == '.py':
            ast.parse(raw.decode('utf-8-sig'), filename=item['path'])
            checked_python += 1
        if path.suffix in {'.py', '.m', '.md', '.json', '.csv', '.txt', '.toml', '.yml', '.tex', '.cff'}:
            text = raw.decode('utf-8-sig')
            if path != Path(__file__).resolve():
                assert not re.search(r'[A-Za-z]:[/\\]+Users[/\\]|/home/[^/\s]+/|/Users/[^/\s]+/', text), item['path']
            assert not re.search(r'gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----', text), item['path']
        if path.suffix == '.json':
            json.loads(raw)
    for item in data_manifest['files']:
        assert hashlib.sha256((ROOT / item['path']).read_bytes()).hexdigest() == item['sha256']
    print(f'PASS: {len(records)} file hashes; {checked_python} Python syntax checks; {len(approved_data)} approved numerical CSV files; no manuscript or generated figures.')
    print('No simulation, numerical test, or plotting job was executed.')


if __name__ == '__main__':
    try:
        main()
    except (AssertionError, ValueError, OSError) as error:
        print(f'FAIL: {error}', file=sys.stderr)
        raise SystemExit(1) from error
