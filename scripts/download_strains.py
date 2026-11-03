"""Acquire the pinned two-strain NCBI experiment inputs."""
import gzip
import json
import os
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from genetrace.io import sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    sources = json.loads((ROOT / 'docs/contig-sources.json').read_text())
    for entry in sources:
        directory = ROOT / 'data/strains' / entry['strain']
        directory.mkdir(parents=True, exist_ok=True)
        for source in entry['files']:
            destination = directory / source['file']
            if not destination.exists():
                partial = destination.with_suffix('.partial')
                with urllib.request.urlopen(source['url'], timeout=60) as response, partial.open('wb') as handle:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
                if sha256(partial) != source['sha256']:
                    raise ValueError('NCBI snapshot changed: ' + source['url'])
                os.replace(partial, destination)
            if sha256(destination) != source['sha256']:
                raise ValueError('Local checksum mismatch: ' + str(destination))
            if str(destination).endswith('.gz'):
                with gzip.open(destination, 'rb') as handle:
                    while handle.read(1024 * 1024):
                        pass
            print(entry['strain'] + '/' + source['file'] + ': verified')


if __name__ == '__main__':
    main()
